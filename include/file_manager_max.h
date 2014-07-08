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

#ifndef __FILE_MANAGER_MAX_H__
#define __FILE_MANAGER_MAX_H__

#include "fasta_file.h"
#include "fastq_file.h"
#include "file_manager.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <string>


class FileManagerMax : public FileManager {
private:
	unsigned long max_nb_reads; // Maximum number of reads per file
public:
	
	// Constructor
	FileManagerMax (const unsigned long & max) {
		current_file = -1;
		total_nb_reads = 0;
		nb_tagged_reads = 0;
		nb_seen_reads = 0;
		max_nb_reads = max;
	}
	
	// Get the next read in the current file or try next file if no more read in current file
	// If the max number of tagged read is reached, then stop
	// If no more file, return an empty read
	std::string & get_next_read_to_filter () {
		std::string & tmp_read = files[current_file]->get_next_read();
		if (nb_tagged_reads >= max_nb_reads) {
			tmp_read.clear();
		}
		while (tmp_read.empty()) {
			current_file++;
			nb_tagged_reads = 0;
			if (current_file >= (int) files.size()) {
				break;
			}
			tmp_read = files[current_file]->get_next_read();
		}
		return tmp_read;
	}
	
	// Get the next read in the current file or try next file if no more read in current file
	// If no more file, return an empty read
	std::string & get_next_read_to_compare () {
		std::string & tmp_read = files[current_file]->get_next_read();
		if (nb_seen_reads >= max_nb_reads) {
			tmp_read.clear();
		}
		while (tmp_read.empty()) {
			current_file++;
			nb_seen_reads = 0;
			if (current_file >= (int) files.size()) {
				break;
			}
			tmp_read = files[current_file]->get_next_read();
		}
		nb_seen_reads++;
		return tmp_read;
	}
	
	void rewind () {
		nb_tagged_reads = 0;
		nb_seen_reads = 0;
		current_file = 0;
		for (std::vector<ReadFile *>::iterator it = files.begin(); it != files.end(); it++) {
			(*it)->rewind();
		}
	}
	
	// Add a file to the FileManager
	void addFile (const std::string & file_name) {
		// Open the given file to check its type (fasta, fastq, gzip ?)
		std::ifstream infile;
		infile.open(file_name.c_str());
		if (!infile.good()) {
			std::cerr << "Cannot open file file " << file_name << " -> ignore\n";
		}
		// Check the first char
		char c = infile.get();
		// Fasta file
		if (c == '>') {
			infile.close();
			files.push_back(new FastaFile(file_name));
		}
		// Fastq file
		else if (c == '@') {
			infile.close();
			files.push_back(new FastqFile(file_name));
		}
		// Gzip file
		else {
			infile.close();
			gzFile tmp_gz_file = (gzFile) gzopen(file_name.c_str(), "r");
			if (!tmp_gz_file) {
				std::cerr << "Cannot open file " << file_name << " -> ignore\n";
				exit(1);
			}
			c = gzgetc(tmp_gz_file);
			// Fasta file
			if (c == '>') {
				gzclose(tmp_gz_file);
				files.push_back(new GzFastaFile(file_name));
			}
			// Fastq file
			else if (c == '@') {
				gzclose(tmp_gz_file);
				files.push_back(new GzFastqFile(file_name));
			} else {
				std::cerr << "Unknown format: " << file_name << " -> ignore\n";
			}
		}
		file_names.push_back (file_name);
		if (sum_of_file_names.empty()) {
			sum_of_file_names = file_name.substr(file_name.rfind("/") + 1);
		} else {
			sum_of_file_names += "-" + file_name.substr(file_name.rfind("/") + 1);
		}
		if (current_file < 0) {
			current_file = 0;
		}
		if (files.back()->nb_valid_reads() < max_nb_reads) {
			total_nb_reads += files.back()->nb_valid_reads();
		} else {
			total_nb_reads += max_nb_reads;
		}
		file_bvs.push_back(files.back()->get_bv());
		file_bvs.back().set_all_false();
	}
	
	// Add a file + boolean vector to the FileManager
	void addFile (const std::string & file_name, const std::string & bv_file_name) {
		// Open the given file to check its type (fasta, fastq, gzip ?)
		std::ifstream infile;
		infile.open(file_name.c_str());
		if (!infile.good()) {
			std::cerr << "Cannot open file " << file_name << " -> ignore\n";
			return;
		}
		// Check the first char
		char c = infile.get();
		if (c == '>') {
			infile.close();
			files.push_back(new FastaFile(file_name, bv_file_name));
		} else if (c == '@') {
			infile.close();
			files.push_back(new FastqFile(file_name, bv_file_name));
		} else {
			infile.close();
			gzFile tmp_gz_file = (gzFile) gzopen(file_name.c_str(), "r");
			if (!tmp_gz_file) {
				std::cerr << "Cannot open file " << file_name << " -> ignore\n";
				return;
			}
			c = gzgetc(tmp_gz_file);
			if (c == '>') {
				gzclose(tmp_gz_file);
				files.push_back(new GzFastaFile(file_name, bv_file_name));
			} else if (c == '@') {
				gzclose(tmp_gz_file);
				files.push_back(new GzFastqFile(file_name, bv_file_name));
			} else {
				std::cerr << "Unknown format: " << file_name << " -> ignore\n";
				return;
			}
		}
		file_names.push_back (file_name);
		sum_of_file_names += file_name.substr(file_name.rfind("/") + 1);
		if (current_file < 0) {
			current_file = 0;
		}
		if (files.back()->nb_valid_reads() < max_nb_reads) {
			total_nb_reads += files.back()->nb_valid_reads();
		} else {
			total_nb_reads += max_nb_reads;
		}
		
		file_bvs.push_back(files.back()->get_bv());
		file_bvs.back().set_all_false();
	}
};

#endif
