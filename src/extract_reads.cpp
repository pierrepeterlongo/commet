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

#include "file_manager.h"

#include <iostream>
#include <string>
#include <ctime>

#include <zlib.h>

std::string version = "2.1";

// -----------------------------------------------------------------------
//                             PROTOTYPES
// -----------------------------------------------------------------------

void print_usage ();

// -----------------------------------------------------------------------
//                                MAIN
// -----------------------------------------------------------------------

int main (int argc, char ** argv)
{
	//const clock_t begin_time = clock();
	
	if (argc<3) {
		print_usage ();
	}
	
	////////////////////////////////////////////////////////////
	// Init parameters
	//
	std::string input_file_name;
	std::string bv_file_name;
	std::string output_file_name;
	bool compress = false;
	
	////////////////////////////////////////////////////////////
	// Read command line arguments
	//
	int arg_pos = 1;
	while (arg_pos < argc){
		std::string flag = argv[arg_pos];
		if (flag[0] != '-') {
			if (input_file_name.empty()) {
				input_file_name = flag;
			} else if (bv_file_name.empty()){
				bv_file_name = flag;
			} else {
				std::cerr << "The mandatory files are already set, unknown file " << flag << " -> ignore\n";
			}
		} else if (flag.compare("-o") == 0) {
			arg_pos++;
			output_file_name = argv[arg_pos];
		} else if (flag.compare("-h") == 0) {
			print_usage ();
			return 0;
		} else if (flag.compare("-v") == 0) {
			std::cout << "\nextract_reads version " << version << "\n";
			return 0;
		} else {
			std::cerr << "Unknown option " << flag << "\n";
            print_usage ();
			return (0);
		}
		arg_pos++;
	}
	
	if (input_file_name.empty()) {
		std::cerr << "Error: An input file name is needed -> exit\n";
		print_usage ();
		return (0);
	} else if (bv_file_name.empty()) {
		std::cerr << "Error: A bv file name is needed -> exit\n";
		print_usage ();
		return (0);
	}
	
	////////////////////////////////////////////////////////////
	// Open the given file to check its type (fasta, fastq, gzip ?)
	//
	ReadFile * read_file = NULL;
	
	std::ifstream infile;
	infile.open(input_file_name.c_str());
	if (!infile.good()) {
		std::cerr << "Cannot open file file " << input_file_name << " -> ignore\n";
	}
	// Check the first char
	std::string basename = input_file_name.substr(input_file_name.rfind("/")+1);
	char c = infile.get();
	if (c == '>') {
		infile.close();
		read_file = new FastaFile(input_file_name, bv_file_name);
	} else if (c == '@') {
		infile.close();
		read_file = new FastqFile(input_file_name, bv_file_name);
	} else {
		infile.close();
		gzFile tmp_gz_file = (gzFile) gzopen(input_file_name.c_str(), "r");
		if (!tmp_gz_file) {
			std::cerr << "Cannot open file " << input_file_name << " -> ignore\n";
			exit(1);
		}
		c = gzgetc(tmp_gz_file);
		if (c == '>') {
			gzclose(tmp_gz_file);
			compress = true;
			read_file = new GzFastaFile(input_file_name, bv_file_name);
		} else if (c == '@') {
			gzclose(tmp_gz_file);
			compress = true;
			read_file = new GzFastqFile(input_file_name, bv_file_name);
		} else {
			std::cerr << "Unknown format: " << input_file_name << " -> ignore\n";
		}
	}
	
	////////////////////////////////////////////////////////////
	// Open the output file and write selected reads in it
	// Different behaviors if : - output file name is given
	//                          - compress is true or not
	//
	if (compress) {
		if (output_file_name.empty()) {
			std::cerr << "Error, try to compress results but no output file name is given\n";
			exit(1);
		}
		gzFile filetmp = (gzFile) gzopen(output_file_name.c_str(), "w6");
		if (filetmp == NULL) {
			std::cerr << "Error, cannot open file " << output_file_name << "\n";
		}
		std::string & current_read = read_file->get_next_read();
		while (!current_read.empty()) {
			gzprintf(filetmp, "%s", read_file->get_data().c_str());
			current_read = read_file->get_next_read();
		}
		gzclose(filetmp);
	} else {
		if (!output_file_name.empty()) {
			std::ofstream outfile;
			outfile.open (output_file_name.c_str());
			if (!outfile.good()) {
				std::cerr << "Cannot write on file " << output_file_name << "\n";
				return 1;
			}
			std::string & current_read = read_file->get_next_read();
			while (!current_read.empty()) {
				outfile << read_file->get_data();
				current_read = read_file->get_next_read();
			}
			outfile.close();
		} else {
			std::string & current_read = read_file->get_next_read();
			while (!current_read.empty()) {
				std::cout << read_file->get_data();
				current_read = read_file->get_next_read();
			}
		}
	}
	
	if (read_file != NULL) {
		delete read_file;
	}
	return 0;
}

void print_usage (){
	std::cout << "\nextract_reads v" << version << "\n";
	std::cout << "Usage:\n\t./extract_reads <input_file> <bv_file> [options]\n";
	std::cout << "Mandatory:\n";
	std::cout << "\t<input_file>\t: file containing reads, in fasta or fastq format, gzipped or not\n";
    std::cout << "\t<bv_file>\t: associated boolean vector file\n";
    std::cout << "Options:\n";
    std::cout << "\t -o string: Output results in the given file [default=stdout]\n";
	std::cout << "\t -h: Prints this message and exit\n";
	std::cout << "\t -v: prints the version number.\n\n";
    exit(0);
}
