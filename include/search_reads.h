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

#ifndef SEARCH_READS_H_
#define SEARCH_READS_H_

#include "bloom_filter.h"
#include "file_manager.h"
#include "alphabet.h"

unsigned long search_reads (const BloomFilter * index, FileManager * search_file_manager, const int & kmer_size, const int & min_hits, unsigned long & nb_searched_reads)
{
	clock_t begin_time = clock();
	// Search reads from search_file_manager in the indexed reads
	HashKey hash (kmer_size);
	Alphabet * alphabet = Alphabet::getInstance();
	
	unsigned long nb_found_reads = 0;
	search_file_manager->rewind();
	std::string & current_read_to_search = search_file_manager->get_next_read_to_compare();
	while (!current_read_to_search.empty()) {
		nb_searched_reads++;
		// Search k-mers
		int seen = 0;
		bool found = false;
		hash.clear();
		for (long i = 0; i < (int) current_read_to_search.size() && !found; i++) {
			if (!alphabet->is_in(current_read_to_search[i])) {
				hash.clear();
			} else if (hash.add(current_read_to_search[i]) >= kmer_size) {
				if (index->is_found(hash)) {
					seen++;
					if (seen >= min_hits) {
						search_file_manager->tag_current_read();
						nb_found_reads++;
						found = true;
					}
					hash.clear();
				}
			}
		}
		// Search the reverse strand if not found in the first step
		if (!found) {
			seen = 0;
			hash.clear();
			for (long i = 0; i < (int) current_read_to_search.size() && !found; i++) {
				if (!alphabet->is_in(current_read_to_search[i])) {
					hash.clear();
				} else if (hash.rv_add(current_read_to_search[i]) >= kmer_size) {
					if (index->is_found(hash)) {
						seen++;
						if (seen >= min_hits) {
							search_file_manager->tag_current_read();
							nb_found_reads++;
							found = true;
						}
						hash.clear();
					}
				}
			}
		}
		current_read_to_search = search_file_manager->get_next_read_to_compare();
	}
	std::cout << float(clock () - begin_time) / CLOCKS_PER_SEC << " seconds\n";
	return nb_found_reads;
}


#endif