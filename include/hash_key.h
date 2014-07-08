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

#ifndef HASH_KEY_H_
#define HASH_KEY_H_

#include <math.h>

class HashKey
{
private:
	unsigned long _keya, _keyb, _keyc, _keyd;
	unsigned long mask_size_kmer;
	unsigned long rv_mask_size_kmer;
	unsigned long bloom_size;
	int hash_size;
public:
	
	HashKey (const int & kmer_size)
	{
		// Set bloom filter size 1000...0
		bloom_size = (unsigned long) pow (2, kmer_size - 1);
		// mask_size_kmer = 1111...1
		mask_size_kmer = (2 * bloom_size) - 1;
		// rv_mask_size_kmer = 0111...1
		rv_mask_size_kmer = bloom_size - 1;
		clear();
	}
	
	// Clear the 4 keys
	void clear()
	{
		hash_size = 0;
		_keya = _keyb = _keyc = _keyd = 0;
	}
	
	// Add an amino acid to the 4 keys and return the hash_size
	//
	// This is supposed to add amino acids from the beginning of the read to the end
	// When an amino acid is added, all bits are pushed to the left
	// then the last bit is set to 1
	const int & add (const char & aa)
	{
		hash_size++;
		_keya = (_keya << 1) & mask_size_kmer;
		_keyb = (_keyb << 1) & mask_size_kmer;
		_keyc = (_keyc << 1) & mask_size_kmer;
		_keyd = (_keyd << 1) & mask_size_kmer;
		if(aa == 'C' || aa == 'c')
		{
			_keyb |= 1;
			_keyc |= 1;
			_keyd |= 1;
		}
		else if(aa == 'G' || aa == 'g')
		{
			_keya |= 1;
			_keyc |= 1;
			_keyd |= 1;
		}
		else if(aa == 'T' || aa == 't')
		{
			_keya |= 1;
			_keyb |= 1;
			_keyd |= 1;
		}
		return hash_size;
	}
	
	// Add the reverse complement amino acid to the 4 keys and return the hash_size
	// This method is designed to be used as the add() method
	// But the produced hash correspond to the reverse complement
	//
	// When an amino acid is added, all bits are pushed to the right
	// then the first bit is set to 1
	const int & rv_add (char & aa)
	{
		hash_size++;
		_keya = (_keya >> 1) & rv_mask_size_kmer;
		_keyb = (_keyb >> 1) & rv_mask_size_kmer;
		_keyc = (_keyc >> 1) & rv_mask_size_kmer;
		_keyd = (_keyd >> 1) & rv_mask_size_kmer;
		if(aa=='A' || aa=='a')
		{
			_keya |= bloom_size;
			_keyb |= bloom_size;
			_keyd |= bloom_size;
		}
		else if(aa == 'C' || aa == 'c')
		{
			_keya |= bloom_size;
			_keyc |= bloom_size;
			_keyd |= bloom_size;
		}
		else if(aa=='G' || aa == 'g')
		{
			_keyb |= bloom_size;
			_keyc |= bloom_size;
			_keyd |= bloom_size;
		}
		return hash_size;
	}
	
	const unsigned long & keya () const {return _keya;}
	const unsigned long & keyb () const {return _keyb;}
	const unsigned long & keyc () const {return _keyc;}
	const unsigned long & keyd () const {return _keyd;}
};

#endif
