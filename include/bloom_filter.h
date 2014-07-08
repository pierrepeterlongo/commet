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

#ifndef BLOOM_FILTER_H_
#define BLOOM_FILTER_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hash_key.h"

/*
 * 4 hash functions
 * keya : A/C -> 0, G/T -> 1
 * keyb : A/G -> 0, C/T -> 1
 * keyc : A/T -> 0, C/G -> 1
 * keyd : A -> 0, C/G/T -> 1
 */

class BloomFilter
{
private:
	char * bloom_vector;
	long bloom_size;
	char MASK_A_EVEN;
	char MASK_B_EVEN;
	char MASK_C_EVEN;
	char MASK_D_EVEN;
	char MASK_A_ODD;
	char MASK_B_ODD;
	char MASK_C_ODD;
	char MASK_D_ODD;

	
public:
	explicit BloomFilter (const int & kmer_size)
	{
		MASK_A_EVEN = 128;
		MASK_B_EVEN = 64;
		MASK_C_EVEN = 32;
		MASK_D_EVEN = 16;
		MASK_A_ODD = 8;
		MASK_B_ODD = 4;
		MASK_C_ODD = 2;
		MASK_D_ODD = 1;
		
		// Set bloom filter size 1000...0
		bloom_size = (unsigned long) pow (2, kmer_size - 1);
		
		// Allocate Bloom filter
		bloom_vector = (char *) calloc (bloom_size, sizeof(char));
		if (bloom_vector == NULL) {
			fprintf(stderr, "Index memory allocation impossible, try with a lower k value or with more RAM memory\n");
			exit(1);
		}
	}
	
	~BloomFilter () {
		free (bloom_vector);
	}
	
	void clear()
	{
		if (bloom_vector != NULL) {
			free (bloom_vector);
		}
		bloom_vector = (char *) calloc (bloom_size, sizeof(char));
		if (bloom_vector == NULL) {
			fprintf(stderr, "Index memory allocation impossible, try with a lower k value or with more RAM memory\n");
			exit(1);
		}
	}
	
	bool empty() {
		for (long i = 0; i < bloom_size; i++) {
			if (bloom_vector[i] > 0) {
				return false;
			}
		}
		return true;
	}
	
	// Add values to the Bloom filter
	// THIS IS WHERE WE PEND MOST OF THE TIME
	// key/2 -> shift the last number, ie : 110010 -> 11001
	// the %2 give if the last number is 0 or 1
	void feed (const HashKey & hash_key)
	{
		bloom_vector[hash_key.keya() / 2] |= (hash_key.keya() % 2 ? MASK_A_ODD : MASK_A_EVEN);
		bloom_vector[hash_key.keyb() / 2] |= (hash_key.keyb() % 2 ? MASK_B_ODD : MASK_B_EVEN);
		bloom_vector[hash_key.keyc() / 2] |= (hash_key.keyc() % 2 ? MASK_C_ODD : MASK_C_EVEN);
		bloom_vector[hash_key.keyd() / 2] |= (hash_key.keyd() % 2 ? MASK_D_ODD : MASK_D_EVEN);
	}

	// Check values to the Bloom filter
	// THIS IS WHERE WE PEND MOST OF THE TIME
	// key/2 -> shift the last number, ie : 110010 -> 11001
	// the %2 give if the last number is 0 or 1
	bool is_found (const HashKey & hash_key) const
	{
		return
		(bloom_vector[hash_key.keya() / 2] & (hash_key.keya() % 2 ? MASK_A_ODD : MASK_A_EVEN)) &&
		(bloom_vector[hash_key.keyb() / 2] & (hash_key.keyb() % 2 ? MASK_B_ODD : MASK_B_EVEN)) &&
		(bloom_vector[hash_key.keyc() / 2] & (hash_key.keyc() % 2 ? MASK_C_ODD : MASK_C_EVEN)) &&
		(bloom_vector[hash_key.keyd() / 2] & (hash_key.keyd() % 2 ? MASK_D_ODD : MASK_D_EVEN));
	}
};

#endif
