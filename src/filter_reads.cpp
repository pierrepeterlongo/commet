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

#include <sstream>
#include <iostream>
#include <cmath>
#include <string>
#include <ctime>
#include <limits.h>
#include "file_manager.h"
#include "alphabet.h"

std::string version = "2.1";

// -----------------------------------------------------------------------
//                             PROTOTYPES
// -----------------------------------------------------------------------

void print_usage ();
int number_of_N (const std::string & read);
float shannon_index(std::string & read);

// -----------------------------------------------------------------------
//                                MAIN
// -----------------------------------------------------------------------

int main (int argc, char ** argv)
{
	const clock_t begin_time = clock();
	
	std::string input_file_name;
	std::string output_file_name;
	int min_size = 0;
    int max_N = INT_MAX;
    float min_shannon = 0.0;
	std::stringstream comment;
	long max_reads = -1;
	long nb_selected_reads = 0;
	
	int arg_pos = 1;
	while (arg_pos < argc){
		std::string flag = argv[arg_pos];
		if (flag[0] != '-') {
			if (input_file_name.empty()) {
				input_file_name = flag;
			} else if (output_file_name.empty()){
				output_file_name = flag;
			} else {
				std::cout << "The mandatory files are already set, unknown file " << flag << " -> ignore\n";
			}
		} else if (flag.compare("-o") == 0) {
			arg_pos++;
			output_file_name = argv[arg_pos];
		} else if (flag.compare("-l") == 0) {
			arg_pos++;
			min_size = atoi (argv[arg_pos]);
		} else if (flag.compare("-n") == 0) {
			arg_pos++;
			max_N = atoi (argv[arg_pos]);
		} else if (flag.compare("-m") == 0) {
			arg_pos++;
			max_reads = atoi (argv[arg_pos]);
		} else if (flag.compare("-e") == 0) {
			arg_pos++;
			min_shannon = atof(argv[arg_pos]);
		} else if (flag.compare("-c") == 0) {
			arg_pos++;
			comment << argv[arg_pos] << "\n";
		} else if (flag.compare("-h") == 0) {
			print_usage ();
			return 0;
		} else if (flag.compare("-v") == 0) {
			std::cout << "\nfilter_reads version " << version << "\n";
			return 0;
		} else {
			std::cerr << "Unknown option " << flag << "\n";
            print_usage ();
			return 1;
		}
		arg_pos++;
	}
	
	if (input_file_name.empty()) {
		std::cerr << "Error: An input file name is needed -> exit\n";
		print_usage ();
		return (0);
	}
	std::string output_message;
	if (output_file_name.empty()) {
		output_message = "No output file name given, results will be written in " + input_file_name + ".bv\n";
		output_file_name = input_file_name + ".bv";
	}
	////////////////////////////////////////////////////////////
	// Open the given file to check its type (fasta, fastq, gzip ?)
	//
	ReadFile * read_file = NULL;
	
	std::ifstream infile;
	infile.open(input_file_name.c_str());
	if (!infile.good()) {
		std::cerr << "Cannot open file " << input_file_name << " -> quit\n";
		return 1;
	}
	// Check the first char
	std::string basename = input_file_name.substr(input_file_name.rfind("/")+1);
	char c = infile.get();
	if (c == '>') {
		infile.close();
		read_file = new FastaFile(input_file_name);
	} else if (c == '@') {
		infile.close();
		read_file = new FastqFile(input_file_name);
	} else {
		infile.close();
		gzFile tmp_gz_file = (gzFile) gzopen(input_file_name.c_str(), "r");
		if (!tmp_gz_file) {
			std::cerr << "Cannot open file " << input_file_name << " -> quit\n";
			exit(1);
		}
		c = gzgetc(tmp_gz_file);
		if (c == '>') {
			gzclose(tmp_gz_file);
			read_file = new GzFastaFile(input_file_name);
		} else if (c == '@') {
			gzclose(tmp_gz_file);
			read_file = new GzFastqFile(input_file_name);
		} else {
			std::cerr << "Unknown format: " << input_file_name << " -> quit\n";
			exit(1);
		}
	}
	
	
	////////////////////////////////////////////////////////////
	// Prepare the comment
	//
	comment << "----------------\n";
	comment << "Reference file\n";
	size_t pos = input_file_name.rfind("/");
	if (pos > 0 && pos < (int) input_file_name.size()) {
		comment << "  " << input_file_name.substr(pos + 1) << "\n";
	} else {
		comment << "  " << input_file_name << "\n";
	}
	comment << "Filter Options\n";
	comment << "  min read size     : " << min_size << "\n";
	if (max_N == INT_MAX) {
		comment << "  max number of N   : infinite\n";
	} else {
		comment << "  max number of N   : " << max_N << "\n";
	}
	
	comment << "  min shannon index : " << min_shannon << "\n";
	
	////////////////////////////////////////////////////////////
	// Test each read to filter values
	//
	if (max_reads == -1) {
		max_reads = read_file->get_nb_reads();
	}
	long nb_rm_length = 0;
	long nb_rm_N = 0;
	long nb_rm_shannon = 0;
	std::string & current_read = read_file->get_next_read();
	while (!current_read.empty() && nb_selected_reads < max_reads) {
		if ((int) current_read.size() < min_size) {
			read_file->untag_current_read();
			nb_rm_length++;
		} else if (number_of_N (current_read) > max_N) {
			read_file->untag_current_read();
			nb_rm_N++;
		} else if (shannon_index(current_read) < min_shannon) {
			read_file->untag_current_read();
			nb_rm_shannon++;
		} else {
			nb_selected_reads++;
		}
		current_read = read_file->get_next_read();
	}
	if (nb_selected_reads >= max_reads) {
		read_file->untag_last_reads();
	}
	read_file->set_bv_comment(comment.str());
	read_file->save_bv(output_file_name);
	
	std::cout << "Length filter [" << min_size << "]: " << nb_rm_length << " reads removed\n";
	if (max_N == INT_MAX) {
		std::cout << "Number of N filter [infinite]: " << nb_rm_N << " reads removed\n";
	} else {
		std::cout << "Number of N filter [" << max_N << "]: " << nb_rm_N << " reads removed\n";
	}
	std::cout << "Shannon filter [" << min_shannon << "]: " << nb_rm_shannon << " reads removed\n";
	std::cout << "Number of selected reads = " << nb_selected_reads << "\n";
	if (!output_message.empty()) {
		std::cout << output_message;
	}
	std::cout << "Total  time : " << float (clock () - begin_time) / CLOCKS_PER_SEC << " s\n";
	return 0;
}


// -----------------------------------------------------------------------
//                                USAGE
// -----------------------------------------------------------------------

void print_usage () {
	std::cout << "\nfilter_reads v" << version << "\n";
	std::cout << "Usage:\n\t./filter_reads <input_file> [options]\n";
	std::cout << "Mandatory:\n";
    std::cout << "\t<input_file>\t: file containing reads, in fasta or fastq format, gzipped or not\n";
    std::cout << "Options:\n";
	std::cout << "\t -o string\t: file where the boolean vector will be written [default=input_file.bv]\n";
    std::cout << "\t -l int\t\t: minimal length a read should have to be kept. [default=0]\n";
    std::cout << "\t -n int\t\t: maximal number of Ns a read should contain to be kept. [default=any]\n";
    std::cout << "\t -e float\t: minimal Shannon index a read should have to be kept. [default=0]\n";
	std::cout << "\t -m int\t\t: maximum number of selected reads [default=all]\n";
    std::cout << "\t -c string\t: the given string will be written in the header of the output file. [default=command line]\n";
	std::cout << "\t -h\t\t: prints this help\n";
	std::cout << "\t -v\t\t: prints the version number.\n\n";
}

// -----------------------------------------------------------------------
//                        NUMBER OF Ns IN A READ
// -----------------------------------------------------------------------

int number_of_N (const std::string & read)
{
	Alphabet * alphabet = Alphabet::getInstance();
	int nb = 0;
	for (int i = 0; i < (int) read.size(); i++) {
		if (!alphabet->is_in(read[i])) {
			nb++;
		}
	}
	return nb;
}

// -----------------------------------------------------------------------
//                        SHANNON INDEX OF A READ
// -----------------------------------------------------------------------

float shannon_index(std::string & read)
{
	float index = 0;
	float freq [5];
	
	// Initiate frequencies to 0.0
	for (int i = 0; i < 5; i++) {
		freq[i] = 0.0;
	}
	
	// Frequency of each letter (A, C, G, T or N)
	for (int i = 0; i < (int) read.size(); i++)
	{
		switch (toupper(read[i]))
		{
			case 'A':
				freq[0] += 1.0;
				break;
			case 'C':
				freq[1] += 1.0;
				break;
			case 'G':
				freq[2] += 1.0;
				break;
			case 'T':
				freq[3] += 1.0;
				break;
			default:
				freq[4] += 1.0;
		}
	}
	
	// Shannon index calculation
	for (int i = 0; i < 5; i++)
	{
		freq[i] = freq[i] / (float) read.size();
		if (freq[i] != 0) {
			index += freq[i] * log (freq[i]) / log(2);
		}
	}
	return fabs (index);
}
