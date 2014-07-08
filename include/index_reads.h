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

#ifndef INDEX_READS_H_
#define INDEX_READS_H_

#include "bloom_filter.h"
#include "file_manager.h"
#include "alphabet.h"

/* index_reads returns a bloom filter of indexed reads from the given file_manager.
 *
 * it also update the number of indexed reads because it may stop if the number of kmer is too high
 *
 * For each read, calculate the hash for each k-mer and feed the bloom filter
 * If an N is found, then reinit hash and continue.
 */
BloomFilter * index_reads (FileManager * index_file_manager, const int & kmer_size, const int & min_hits, const unsigned long & max_kmer, unsigned long & nb_indexed_reads)
{
	unsigned long nb_indexed_kmers = 0;
	BloomFilter * bloom_filter = new BloomFilter (kmer_size);
	HashKey hash (kmer_size);
	Alphabet * alphabet = Alphabet::getInstance();
	
	std::string & current_read_to_index = index_file_manager->get_next_read_to_compare();
	while (!current_read_to_index.empty() && nb_indexed_kmers < max_kmer) {
		nb_indexed_reads++;
		hash.clear();
		for (int i = 0; i < (int) current_read_to_index.size(); i++) {
			if (!alphabet->is_in(current_read_to_index[i])) {
				hash.clear();
			} else if (hash.add(current_read_to_index[i]) >= kmer_size) {
				bloom_filter->feed(hash);
				nb_indexed_kmers++;
			}
		}
		current_read_to_index = index_file_manager->get_next_read_to_compare();
	}
	return bloom_filter;
}

#endif
