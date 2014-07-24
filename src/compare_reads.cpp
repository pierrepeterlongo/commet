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
#include "set_parser.h"


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

// -----------------------------------------------------------------------
//                                MAIN
// -----------------------------------------------------------------------

int main (int argc, char ** argv)
{
	// A files variables
	std::string A_file_list;
	FileManager * A_set = new FileManager;
	std::map <std::string, std::vector <std::string> > A_file_names;
	std::map <std::string, std::vector <std::string> > A_bv_names;
	
	// B files variables
	std::string B_file_list;
	FileManager * B_set = new FileManager;
	std::map <std::string, std::vector <std::string> > B_file_names;
	std::map <std::string, std::vector <std::string> > B_bv_names;
	
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
			// File name for list of A files
			arg_pos++;
			if (arg_pos >= argc) {
				std::cerr << "Error, flag " << argv[arg_pos - 1] << " needs an argument\n";
				print_usage();
				exit(1);
			}
			if (!A_file_list.empty()) {
				std::cerr << "A files already given (-i) -> ignore";
			} else {
				A_file_list = argv[arg_pos];
			}
		} else if (flag.compare("-s") == 0) {
			// File name for list of B files
			arg_pos++;
			if (arg_pos >= argc) {
				std::cerr << "Error, flag " << argv[arg_pos - 1] << " needs an argument\n";
				print_usage();
				exit(1);
			}
			if (!B_file_list.empty()) {
				std::cerr << "B files already given (-s) -> ignore";
			} else {
				B_file_list = argv[arg_pos];
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
			std::cout << "\ncompare_reads version " << version << "\n";
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
	// Put A files in a file manager
	//
	read_sets(A_file_list, A_file_names, A_bv_names);
	if (A_file_names.size() != 1) {
		std::cerr << "Only one set of files is allowed for A -> keep first set only\n";
	}
	A_set->set_nickname(A_file_names.begin()->first);
	std::vector <std::string> tmp_A_file_names = A_file_names.begin()->second;
	std::vector <std::string> tmp_A_bv_names = A_bv_names.begin()->second;
	for (size_t file_pos = 0; file_pos < tmp_A_file_names.size(); file_pos++) {
		if (tmp_A_bv_names[file_pos].empty()) {
			std::cout << "open " << tmp_A_file_names[file_pos] << "\n";
			A_set->addFile(tmp_A_file_names[file_pos]);
		} else {
			std::cout << "open " << tmp_A_file_names[file_pos] << "," << tmp_A_bv_names[file_pos] << "\n";
			A_set->addFile(tmp_A_file_names[file_pos], tmp_A_bv_names[file_pos]);
		}
	}
	
	////////////////////////////////////////////////////////////
	// Put B files in B file manager
	//
	read_sets(B_file_list, B_file_names, B_bv_names);
	if (B_file_names.size() != 1) {
		std::cerr << "Only one set of files is allowed for B -> keep first set only\n";
	}
	B_set->set_nickname(B_file_names.begin()->first);
	std::vector <std::string> tmp_B_file_names = B_file_names.begin()->second;
	std::vector <std::string> tmp_B_bv_names = B_bv_names.begin()->second;
	for (size_t file_pos = 0; file_pos < tmp_B_file_names.size(); file_pos++) {
		if (tmp_B_bv_names[file_pos].empty()) {
			std::cout << "open " << tmp_B_file_names[file_pos] << "\n";
			B_set->addFile(tmp_B_file_names[file_pos]);
		} else {
			std::cout << "open " << tmp_B_file_names[file_pos] << "," << tmp_B_bv_names[file_pos] << "\n";
			B_set->addFile(tmp_B_file_names[file_pos], tmp_B_bv_names[file_pos]);
		}
	}
	
	////////////////////////////////////////////////////////////
	// Create the index in a BloomFilter
	// and
	// Search files
	//
	BloomFilter * index = NULL;
	unsigned long nb_reads_A = A_set->get_total_nb_reads();
	unsigned long nb_reads_B = B_set->get_total_nb_reads();
	
	////////////////////////////////////////////////////////////
	// search B in A
	//
	std::cout << "\n------------------------------------------------------------------\n";
	std::cout << "finding reads from {" << B_file_names.begin()->first << "} present in raw {" << A_file_names.begin()->first << "}\n";
	std::cout << "------------------------------------------------------------------\n";
	unsigned long nb_reads_to_index = A_set->get_total_nb_reads();
	unsigned long nb_indexed_reads = 0;
	unsigned long nb_found_reads = 0;
	unsigned long nb_searched_reads = 0;
	clock_t index_time = 0;
	clock_t search_time = 0;
	clock_t start_time = clock();
	while (nb_indexed_reads < nb_reads_to_index) {
		if (index != NULL) {
			delete index;
			index = NULL;
		}
		const clock_t index_start = clock();
		index = index_reads (A_set, kmer_size, min_hits, max_kmer, nb_indexed_reads);
		index_time += clock() - index_start;
		const clock_t search_start = clock();
		nb_found_reads += search_reads(index, B_set, kmer_size, min_hits, nb_searched_reads);
		search_time += clock() - search_start;
	}
	B_set->apply_bv_on_files();
	std::cout << "Index  time: " << float (index_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "Search time: " << float (search_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "Total  time: " << float (clock() - start_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "[indexed " << nb_indexed_reads << ", searched " << nb_searched_reads << ", shared " << nb_found_reads << "]\n";
	
	////////////////////////////////////////////////////////////
	// search A in (B in A)
	//
	std::cout << "\n------------------------------------------------------------------\n";
	std::cout << "finding reads from {" << A_file_names.begin()->first << "} present in raw {" << B_file_names.begin()->first << "} present in raw {" << A_file_names.begin()->first << "}\n";
	std::cout << "------------------------------------------------------------------\n";
	nb_reads_to_index = B_set->get_total_nb_reads();
	nb_indexed_reads = 0;
	nb_found_reads = 0;
	nb_searched_reads = 0;
	index_time = 0;
	search_time = 0;
	B_set->rewind();
	A_set->rewind();
	start_time = clock();
	while (nb_indexed_reads < nb_reads_to_index) {
		if (index != NULL) {
			delete index;
			index = NULL;
		}
		const clock_t index_start = clock();
		index = index_reads (B_set, kmer_size, min_hits, max_kmer, nb_indexed_reads);
		index_time += clock() - index_start;
		const clock_t search_start = clock();
		nb_found_reads += search_reads(index, A_set, kmer_size, min_hits, nb_searched_reads);
		search_time += clock() - search_start;
	}
	A_set->save_bv(out_path, B_set->get_nickname());
	A_set->apply_bv_on_files();
	std::cout << "Index  time: " << float (index_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "Search time: " << float (search_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "Total  time: " << float (clock() - start_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "[indexed " << nb_indexed_reads << ", searched " << nb_searched_reads << ", shared " << nb_found_reads << "] "<< 100 * (float) nb_found_reads / (float) nb_reads_A <<"%\n";
	
	////////////////////////////////////////////////////////////
	// search B in (A in (B in A))
	//
	std::cout << "\n------------------------------------------------------------------\n";
	std::cout << "finding reads from {" << B_file_names.begin()->first << "} present in raw {" << A_file_names.begin()->first << "} present in raw {" << B_file_names.begin()->first << "} present in raw {" << A_file_names.begin()->first << "}\n";
	std::cout << "------------------------------------------------------------------\n";
	nb_reads_to_index = A_set->get_total_nb_reads();
	nb_indexed_reads = 0;
	nb_found_reads = 0;
	nb_searched_reads = 0;
	index_time = 0;
	search_time = 0;
	start_time = clock();
	B_set->rewind();
	A_set->rewind();
	while (nb_indexed_reads < nb_reads_to_index) {
		if (index != NULL) {
			delete index;
			index = NULL;
		}
		const clock_t index_start = clock();
		index = index_reads (A_set, kmer_size, min_hits, max_kmer, nb_indexed_reads);
		index_time += clock() - index_start;
		const clock_t search_start = clock();
		nb_found_reads += search_reads(index, B_set, kmer_size, min_hits, nb_searched_reads);
		search_time += clock() - search_start;
	}
	B_set->save_bv(out_path, A_set->get_nickname());
	std::cout << "Index  time: " << float (index_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "Search time: " << float (search_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "Total  time: " << float (clock() - start_time) / CLOCKS_PER_SEC << " s\n";
	std::cout << "[indexed " << nb_indexed_reads << ", searched " << nb_searched_reads << ", shared " << nb_found_reads << "] "<< 100 * (float) nb_found_reads / (float) nb_reads_B <<"%\n";
	if (index != NULL) {
		delete index;
		index = NULL;
	}
	return 0;
}

// -----------------------------------------------------------------------
//                             PRINT USAGE
// -----------------------------------------------------------------------
void print_usage ()
{
	std::cerr << "\ncompare_reads, version " << version << "\n";
	std::cerr << "Usage : ./compare_reads -i <file> -s <file> [options]\n";
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
