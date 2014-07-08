/*
 * Contributors :
 *   Pierre PETERLONGO, pierre.peterlongo@inria.fr [12/06/13]
 *   Nicolas MAILLET, nicolas.maillet@inria.fr     [12/06/13]
 *   Guillaume Collet, guillaume@gcollet.fr        [27/05/14]
 *
 * This software is a computer program whose purpose is to find all the
 * similar reads between two set of NGS reads. It also provide a similarity
 * score between the two samples.
 *
 * Copyright (C) 2014  INRIA
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "search_reads.h" 
#include "index_reads.h"
#include "file_manager.h"
#include "file_manager_max.h"
#include "bloom_filter.h"
#include "boolean_vector.h"


#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <map>

std::string version = "2.1";

// -----------------------------------------------------------------------
//                              PROTOTYPE
// -----------------------------------------------------------------------
void print_usage ();

void read_sets (const std::string & file_name, std::map <std::string, std::vector <std::string> > & file_names, std::map <std::string, std::vector <std::string> > & bv_names);
void remove_spaces (std::string & fname);

// -----------------------------------------------------------------------
//                                MAIN
// -----------------------------------------------------------------------

int main (int argc, char ** argv)
{
	// search files variables
	std::string search_file_list;
	std::vector <FileManager * > search_sets;
	std::map <std::string, std::vector <std::string> > search_file_names;
	std::map <std::string, std::vector <std::string> > search_bv_names;
	
	// index files variables
	std::string index_file_list;
	FileManager * index_set = new FileManager;
	std::map <std::string, std::vector <std::string> > index_file_names;
	std::map <std::string, std::vector <std::string> > index_bv_names;
	
	// general parameters
	int kmer_size = 33;
	int min_hits = 2;
	unsigned long max_kmer = (unsigned long) (1000000000.0 / pow (2, 33 - kmer_size));
	
	// path to output messages
	std::string log_path = ".";
	std::string out_path = ".";
	
	////////////////////////////////////////////////////////////
	// Read command line arguments
	//
	if (argc == 1) {
		print_usage ();
		return (0);
	}
	std::string tmp;
	int arg_pos = 1;
	while (arg_pos < argc){
		std::string flag = argv[arg_pos];
		if (flag.compare("-i") == 0) {
			// File name for list of index files
			arg_pos++;
			if (arg_pos >= argc) {
				std::cerr << "Error, flag " << argv[arg_pos - 1] << " needs an argument\n";
				print_usage();
				exit(1);
			}
			if (!index_file_list.empty()) {
				std::cerr << "index files already given (-i) -> ignore";
			} else {
				index_file_list = argv[arg_pos];
			}
		} else if (flag.compare("-s") == 0) {
			// File name for list of search files
			arg_pos++;
			if (arg_pos >= argc) {
				std::cerr << "Error, flag " << argv[arg_pos - 1] << " needs an argument\n";
				print_usage();
				exit(1);
			}
			if (!search_file_list.empty()) {
				std::cerr << "search files already given (-s) -> ignore";
			} else {
				search_file_list = argv[arg_pos];
			}
		} else if (flag.compare("-l") == 0) {
			// The path to the log directory
			arg_pos++;
			if (arg_pos >= argc) {
				std::cerr << "Error, flag " << argv[arg_pos - 1] << " needs an argument\n";
				print_usage();
				exit(1);
			}
			log_path = argv[arg_pos];
		} else if (flag.compare("-o") == 0) {
			// The path to the ouput directory
			arg_pos++;
			if (arg_pos >= argc) {
				std::cerr << "Error, flag " << argv[arg_pos - 1] << " needs an argument\n";
				print_usage();
				exit(1);
			}
			out_path = argv[arg_pos];
		} else if (flag.compare("-k") == 0) {
			// The size of k-mers
			arg_pos++;
			if (arg_pos >= argc) {
				std::cerr << "Error, flag " << argv[arg_pos - 1] << " needs an argument\n";
				print_usage();
				exit(1);
			}
			kmer_size = atoi(argv[arg_pos]);
			max_kmer = (unsigned long) (1000000000.0 / pow (2, 33 - kmer_size));
			std::cout << "k-mer size (-k) = " << kmer_size << "\n";
		} else if (flag.compare("-t") == 0) {
			// The minimal number of hits to consider two reads as similar
			arg_pos++;
			if (arg_pos >= argc) {
				std::cerr << "Error, flag " << argv[arg_pos - 1] << " needs an argument\n";
				print_usage();
				exit(1);
			}
			min_hits = atoi(argv[arg_pos]);
			std::cout << "min hits (-t) = " << min_hits << "\n";
		} else if (flag.compare("-h") == 0) {
			print_usage ();
			return 0;
		} else if (flag.compare("-v") == 0) {
			std::cout << "\nindex_and_search version " << version << "\n";
			return 0;
		} else {
			std::cerr << "Unknown option " << flag << "\n";
			print_usage ();
			return (0);
		}
		arg_pos++;
	}
	
	////////////////////////////////////////////////////////////
	// Check existence of out_path and log_path
	// if not then create them
	//
	struct stat info;
	if (stat (log_path.c_str(), &info ) != 0){
		mkdir(log_path.c_str(), S_IRWXU|S_IRGRP|S_IXGRP);
	} else if (!(info.st_mode & S_IFDIR)){
		std::cerr << "Error: " << log_path << " already exists and is not a directory\n";
		exit(1);
	}
	
	if (stat (out_path.c_str(), &info ) != 0){
		mkdir(out_path.c_str(), S_IRWXU|S_IRGRP|S_IXGRP);
	} else if (!(info.st_mode & S_IFDIR)) {
		std::cerr << "Error: " << out_path << " already exists and is not a directory\n";
		exit(1);
	}
	
	////////////////////////////////////////////////////////////
	// Put index files in a file manager
	//
	read_sets(index_file_list, index_file_names, index_bv_names);
	if (index_file_names.size() != 1) {
		std::cerr << "Only one set of files is allowed for indexing\n";
		exit(1);
	}
	index_set->set_nickname(index_file_names.begin()->first);
	std::vector <std::string> tmp_index_file_names = index_file_names.begin()->second;
	std::vector <std::string> tmp_index_bv_names = index_bv_names.begin()->second;
	for (size_t file_pos = 0; file_pos < tmp_index_file_names.size(); file_pos++) {
		if (tmp_index_bv_names[file_pos].empty()) {
			std::cout << "open " << tmp_index_file_names[file_pos] << "\n";
			index_set->addFile(tmp_index_file_names[file_pos]);
		} else {
			std::cout << "open " << tmp_index_file_names[file_pos] << "," << tmp_index_bv_names[file_pos] << "\n";
			index_set->addFile(tmp_index_file_names[file_pos], tmp_index_bv_names[file_pos]);
		}
	}
	
	////////////////////////////////////////////////////////////
	// Put search files in search file manager
	//
	read_sets (search_file_list, search_file_names, search_bv_names);
	for (std::map <std::string, std::vector <std::string> >::iterator it_set = search_file_names.begin(); it_set != search_file_names.end(); it_set++) {
		FileManager * current_manager = new FileManager ();
		current_manager->set_nickname(it_set->first);
		for (size_t file_pos = 0; file_pos < it_set->second.size(); file_pos++) {
			if (search_bv_names[it_set->first][file_pos].empty()) {
				std::cout << "open " << it_set->second[file_pos] << "\n";
				current_manager->addFile(it_set->second[file_pos]);
			} else {
				std::cout << "open " << it_set->second[file_pos] << "," << search_bv_names[it_set->first][file_pos] << "\n";
				current_manager->addFile(it_set->second[file_pos],search_bv_names[it_set->first][file_pos]);
			}
		}
		search_sets.push_back(current_manager);
	}
	
	////////////////////////////////////////////////////////////
	// Create the index in a BloomFilter
	// and
	// Search files in search_sets
	//
	unsigned long nb_reads_to_index = index_set->get_total_nb_reads();
	unsigned long nb_indexed_reads = 0;
	std::vector<unsigned long> nb_found_reads (search_sets.size(), 0);
	std::vector<unsigned long> nb_searched_reads (search_sets.size(), 0);
	
	BloomFilter * index = NULL;
	
	while (index_set->get_reads_count() < nb_reads_to_index) {
		if (index != NULL) {
			delete index;
			index = NULL;
		}
		// index
		index = index_reads (index_set, kmer_size, min_hits, max_kmer, nb_indexed_reads);
		
		// search
		for (size_t set_pos = 0; set_pos < search_sets.size(); set_pos++) {
			nb_found_reads[set_pos] += search_reads(index, search_sets[set_pos], kmer_size, min_hits, nb_searched_reads[set_pos]);
			std::cout << search_sets[set_pos]->get_nickname() << " in " << index_set->get_nickname() << "\n" << "[indexed " << nb_indexed_reads << ", searched " << nb_searched_reads[set_pos] << ", shared " << nb_found_reads[set_pos] << "]\n";
		}
	}
	if (index != NULL) {
		delete index;
	}
	
	for (size_t set_pos = 0; set_pos < search_sets.size(); set_pos++) {
		search_sets[set_pos]->save_bv(out_path, index_file_names.begin()->first);
	}
	return 0;
}

// -----------------------------------------------------------------------
//                             PRINT USAGE
// -----------------------------------------------------------------------
void print_usage ()
{
	std::cerr << "\nindex_and_search, version " << version << "\n";
	std::cerr << "Usage : ./index_and_search -i <file> -s <file> [options]\n";
	std::cerr << "Mandatory:\n";
	std::cerr << "\t -i <file>: A file containing the list of files to index (comma separated) - MANDATORY\n";
	std::cerr << "\t            Each line of the file corresponds to a set of files (comma separated)\n";
	std::cerr << "\t -s <file>: A file containing the list of file sets to search - MANDATORY\n";
	std::cerr << "\t            Each line of the file corresponds to a set of files (comma separated)\n";
	std::cerr << "Options:\n";
	std::cerr << "\t -l </.../>: ABSOLUTE path to log folder\n";
	std::cerr << "\t -o </.../>: ABSOLUTE path to output folder\n";
	std::cerr << "\t -k <value>: Size of k-mers (value of k). [default=32]\n";
	std::cerr << "\t -t <value>: Number of shared k-mers. [default=2]\n";
	std::cerr << "\t -h: Prints this message and exit\n";
	std::cerr << "\t -v: Prints the version number and exit\n";
}

void remove_spaces (std::string & fname)
{
	while (fname[0] == ' ') {
		fname = fname.substr(1);
	}
	while (fname[fname.size()-1] == ' ') {
		fname = fname.substr(0, fname.size() - 1);
	}
}

// -----------------------------------------------------------------------
//                               READ_SETS
// -----------------------------------------------------------------------
void read_sets (const std::string & file_name, std::map <std::string, std::vector <std::string> > & file_names, std::map <std::string, std::vector <std::string> > & bv_names)
{
	file_names.clear();
	bv_names.clear();
	
	std::ifstream infile;
	infile.open(file_name.c_str());
	if (!infile.good()) {
		std::cerr << "Cannot read file " << file_name << "\n";
		exit(1);
	}
	int nb_sets = 0;
	while (infile.good()) {
		std::string line;
		getline(infile, line);
		if (!line.empty()) {
			nb_sets++;
			std::stringstream current_tag;
			if (line.find(":") < line.size()) {
				current_tag << line.substr(0, line.find(":"));
				line = line.substr(line.find(":") + 1);
			} else {
				current_tag << "SET" << nb_sets;
			}
			std::vector <std::string> current_files;
			std::vector <std::string> current_bvs;
			while (!line.empty() && line.find(";") < line.size()) {
				std::string fname = line.substr(0, line.find(";"));
				remove_spaces(fname);
				std::string bv_name;
				if (fname.find(",") < fname.size()) {
					bv_name = fname.substr(fname.find(",") + 1);
					remove_spaces(bv_name);
					fname = fname.substr(0,fname.find(","));
					remove_spaces(fname);
				}
				current_files.push_back(fname);
				current_bvs.push_back(bv_name);
				line = line.substr(line.find(";") + 1);
			}
			std::string fname = line;
			remove_spaces(fname);
			std::string bv_name;
			if (fname.find(",") < fname.size()) {
				bv_name = fname.substr(fname.find(",") + 1);
				remove_spaces(bv_name);
				fname = fname.substr(0,fname.find(","));
				remove_spaces(fname);
			}
			current_files.push_back(fname);
			current_bvs.push_back(bv_name);
			file_names[current_tag.str()] = current_files;
			bv_names[current_tag.str()] = current_bvs;
		}
	}
	infile.close();
}

