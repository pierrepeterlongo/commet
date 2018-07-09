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

#ifndef BOOLEAN_VECTOR_H_
#define BOOLEAN_VECTOR_H_

#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <fcntl.h>
#include <string>

/*
 * The BooleanVector class implements an array of bits
 */
class BooleanVector
{
private:
	// contant char = 11111111
	static const char full_char = (char) 255;
	
	// the array of bits is an array of char (8 bits per char)
	char * boolean_vector;
	
	// mask to access the 8 different bits of a char independantly
	char mask [8];
	
	// size of the boolean vector (in bits)
	unsigned long boolean_vector_size;
	
	// size of the boolean vector (in char)
	unsigned long boolean_vector_char_size;
	
	// The comment at the beginning of the file
	std::string comment;
	
public:
	
	//
	// Construct an empty boolean vector
	//
	BooleanVector ()
	{
		mask[0]=1;   //00000001
		mask[1]=2;   //00000010
		mask[2]=4;   //00000100
		mask[3]=8;   //00001000
		mask[4]=16;  //00010000
		mask[5]=32;  //00100000
		mask[6]=64;  //01000000
		mask[7]=128; //10000000
		boolean_vector = NULL;
		boolean_vector_size = 0;
		boolean_vector_char_size = 0;
	};
	
	//
	// Construct by copy
	//
	BooleanVector (const BooleanVector & bv)
	{
		mask[0]=1;   //00000001
		mask[1]=2;   //00000010
		mask[2]=4;   //00000100
		mask[3]=8;   //00001000
		mask[4]=16;  //00010000
		mask[5]=32;  //00100000
		mask[6]=64;  //01000000
		mask[7]=128; //10000000
		boolean_vector = NULL;
		boolean_vector_size = bv.size();
		boolean_vector_char_size = boolean_vector_size / 8 + 1;
		boolean_vector = (char *) calloc (boolean_vector_char_size, sizeof(char));
		if (boolean_vector == NULL) {
			std::cerr << "Cannot allocate memory for variable boolean_vector, exit\n";
			exit(1);
		}
		memcpy(boolean_vector, bv.get_vector(), boolean_vector_char_size);
	}
	
	//
	// Destruct the boolean vector
	//
	~BooleanVector ()
	{
		free(boolean_vector);
	}
	
	//
	// Copy operator
	//
	BooleanVector operator = (const BooleanVector & bv)
	{
		init_false(bv.size());
		memcpy(boolean_vector, bv.get_vector(), boolean_vector_char_size);
		return * this;
	}
	
	//
	// Initiate the boolean vector of the given size (all bits are 0)
	//
	void init_false (const unsigned long & size) {
		boolean_vector_size = size;
		boolean_vector_char_size = boolean_vector_size / 8 + 1;
		if (boolean_vector != NULL) {
			free(boolean_vector);
			boolean_vector = NULL;
		}
		boolean_vector = (char *) calloc (boolean_vector_char_size, sizeof(char));
		if (boolean_vector == NULL) {
			std::cerr << "Cannot allocate memory for variable boolean_vector, exit\n";
			exit(1);
		}
	}
	
	//
	// Initiate the boolean vector of the given size (all bits are 1)
	//
	void init_true (const unsigned long & size) {
		boolean_vector_size = size;
		boolean_vector_char_size = boolean_vector_size / 8 + 1;
		if (boolean_vector != NULL) {
			free(boolean_vector);
			boolean_vector = NULL;
		}
		boolean_vector = (char *) malloc (boolean_vector_char_size * sizeof(char));
		if (boolean_vector == NULL) {
			std::cerr << "Cannot allocate memory for variable boolean_vector, exit\n";
			exit(1);
		}
		memset (boolean_vector, 255, boolean_vector_char_size);
		for (size_t i = boolean_vector_size; i < boolean_vector_char_size * 8; i++) {
			unset(i);
		}
	}
    
    // Generates a bv of a given length with a given proportion of 1's
    void random_vector(const float percentage) {
        /* initialize random seed: */
        srand (time(NULL));
        init_false(boolean_vector_size); 
        for (unsigned long i=0; i<boolean_vector_size;i++){
            if (rand() % 100000 < 1000*percentage) set(i);
        }
    }
	
	//
	// Remove the boolean vector and set sizes to 0
	//
	void clear ()
	{
		free(boolean_vector);
		boolean_vector = NULL;
		boolean_vector_size = 0;
		boolean_vector_char_size = 0;
	}
	
	//
	// Reinit the boolean vector to 0
	//
	void set_all_false ()
	{
		if (boolean_vector != NULL) {
			free(boolean_vector);
			boolean_vector = NULL;
		}
		boolean_vector = (char *) calloc (boolean_vector_char_size, sizeof(char));
		if (boolean_vector == NULL) {
			std::cerr << "Cannot allocate memory for variable boolean_vector, exit\n";
			exit(1);
		}
	}
	
	//
	// Reinit the boolean vector to 1
	//
	void set_all_true ()
	{
		memset (boolean_vector, 255, boolean_vector_char_size);
	}

	
	//
	// Test the size of the boolean vector
	//
	bool empty () {
		return boolean_vector_size == 0;
	}
	
	//
	// Test if the bit at position i is 1
	//
	char is_set (const unsigned long & i) const
	{
		return (boolean_vector[i / 8] & mask[i % 8]);
	}
	
	//
	// Set the bit at position i to 1
	//
	void set (const unsigned long & i)
	{
		boolean_vector[i / 8] |= mask[i % 8];
	}
	
	//
	// Set the bit at position i to 1
	//
	void unset (const unsigned long & i)
	{
		boolean_vector[i / 8] &= (~ mask[i % 8]);
	}
	
	// Get the number of bits to 1
	unsigned long nb_one ()
	{
		unsigned long res = 0;
		unsigned long i = 0;
		for (; i < boolean_vector_char_size - 1; i++) {
			if (boolean_vector[i]) {
				if (boolean_vector[i] == full_char) {
					res += 8;
				} else {
					char tmp = boolean_vector[i];
					tmp = (0x55 & tmp) + (0x55 & (tmp >> 1));
					tmp = (0x33 & tmp) + (0x33 & (tmp >> 2));
					tmp = (0x0f & tmp) + (0x0f & (tmp >> 4));
					res += (int) tmp;
				}
			}
		}
		char tmp = boolean_vector[boolean_vector_char_size - 1];
		tmp = (0x55 & tmp) + (0x55 & (tmp >> 1));
		tmp = (0x33 & tmp) + (0x33 & (tmp >> 2));
		tmp = (0x0f & tmp) + (0x0f & (tmp >> 4));
		res += (int) tmp;
		if (res > boolean_vector_size) {
			res = boolean_vector_size;
		}
		return res;
	}

		
	//
	// Get the size of the boolean vector
	//
	const unsigned long & size () const
	{
		return boolean_vector_size;
	}
	
	// Get the boolean vector
	const char * get_vector () const {
		return boolean_vector;
	}
	
	// Write the boolean vector on stdout in human readable format
	void print () const
	{
		// Print comment + size
		std::cout << comment << "\n#" << boolean_vector_size << "\n";
		// Print the boolean values
		for (unsigned long int i = 0; i < boolean_vector_char_size; i++) {
			std::cout << boolean_vector[i];
		}
	}
	
	//
	// Write the boolean vector in the given file
	// The first bytes store the boolean_vector_size (unsigned long)
	// then the boolean_vector array of char.
	//
	void print (const std::string & file_name) const
	{
		// Prepare comment + size
		std::stringstream tmp_str;
		tmp_str << comment << "\n#" << boolean_vector_size << "\n";
		
		// Open file for writing
		int fd = open (file_name.c_str(), O_RDWR | O_CREAT | O_TRUNC, (mode_t) 0600);
		if (fd == -1) {
			std::cerr << "Error opening file " << file_name << " -> exit\n";
			exit(1);
		}
		// Resize the file to fit the boolean vector
		unsigned long file_size = tmp_str.str().size() * sizeof(char) + boolean_vector_char_size;
		long long result = lseek (fd, file_size - 1, SEEK_SET);
		if (result == -1) {
			std::cerr << "Error resizing file " << file_name << " -> exit\n";
			close(fd);
			exit(1);
		}
		// Write a dumb byte at the end of the file
		result = write (fd, "", 1);
		if (result != 1) {
			std::cerr << "Error writing last byte of " << file_name << " -> exit\n";
			close(fd);
			exit(1);
		}
		// map the file
		char * map = (char *) mmap(0, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
		if (map == MAP_FAILED) {
			std::cerr << "Error mapping file " << file_name << " -> exit\n";
			close(fd);
			exit(1);
		}
		// copy the comment + size part
		memcpy (map, tmp_str.str().c_str(), tmp_str.str().size());
		// copy the boolean vector directly in the file
		memcpy (map + tmp_str.str().size(), boolean_vector, boolean_vector_char_size);
		// unmap the file
		if (munmap(map, file_size) == -1) {
			std::cerr << "Error un-mapping file " << file_name << " -> exit\n";
			close (fd);
		}
		close(fd);
	}
	
	//
	// Read a boolean vector in the given file
	// read the first bytes to get the size of the boolean vector
	// then allocate the memory and copy the values of the array
	//
	void read (const std::string & file_name)
	{
		// open the file
		int fd = open (file_name.c_str(), O_RDONLY);
		if (fd == -1) {
			std::cerr << "Error opening file " << file_name << " -> exit\n";
			exit(1);
		}
		// get statistics to find the size of the file
		struct stat sb;
		if (fstat (fd, &sb) == -1) {
			std::cerr << "Error getting statistics from file " << file_name << " -> exit\n";
			close(fd);
			exit(1);
		}
		// map the file
		char * map = (char *) mmap (0, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
		if (map == MAP_FAILED) {
			std::cerr << "Error mapping file " << file_name << " -> exit\n";
			close(fd);
			exit(1);
		}
		// If the boolean vector was already allocated, erase data
		if (boolean_vector != NULL) {
			free(boolean_vector);
			boolean_vector = NULL;
		}
		if (!comment.empty()) {
			comment.clear();
		}
		long i = 0;
		while (map[i] != '#' && i < sb.st_size) {
			comment += map[i];
			i++;
		}
		i++;
		comment = comment.substr(0, comment.size() - 1);
		std::string tmp_str;
		while (map[i] != '\n' && i < sb.st_size) {
			tmp_str += map[i];
			i++;
		}
		i++;
		if (tmp_str.empty()) {
			std::cerr << "Error, boolean vector does not contain its size\n";
			exit(1);
		}
		unsigned long size = atoi(tmp_str.c_str());
		// Allocate memory for the boolean vector
		init_false (size);
		
		// Directly copy the boolean vector
		memcpy (boolean_vector, &map[i], boolean_vector_char_size);
		// unmap the file
		if (munmap (map, sb.st_size) == -1) {
			std::cerr << "Error un-mapping file " << file_name << " -> exit\n";
			close(fd);
			exit(1);
		}
		// close the file
		close(fd);
	}
	
		
	// Apply a 'and' operator between this boolean vector and bv2
	void full_and (BooleanVector & bv2)
	{
		if (bv2.size() != boolean_vector_size) {
			std::cerr << "Error: the two vectors are not the same size -> exit\n";
			exit(1);
		}
		const char * bv2_vector = bv2.get_vector();
		for (unsigned long i = 0; i < boolean_vector_char_size; i++) {
			boolean_vector[i] = boolean_vector[i] & bv2_vector[i];
		}
	}
	
	// Apply a 'or' operator between this boolean vector and bv2
	void full_or (BooleanVector & bv2)
	{
		if (bv2.size() != boolean_vector_size) {
			std::cerr << "Error: the two vectors are not the same size -> exit\n";
			exit(1);
		}
		const char * bv2_vector = bv2.get_vector();
		for (unsigned long i = 0; i < boolean_vector_char_size; i++) {
			boolean_vector[i] = boolean_vector[i] | bv2_vector[i];
		}
	}
	
	// Apply a 'not' operator on this boolean vector
	void full_not ()
	{
		for (unsigned long i = 0; i < boolean_vector_char_size; i++) {
			boolean_vector[i] = ~boolean_vector[i];
		}
	}
	
	// Apply a 'and not' operator between this boolean vector and bv2
	void full_and_not (BooleanVector & bv2)
	{
		if (bv2.size() != boolean_vector_size) {
			std::cerr << "Error: the two vectors are not the same size -> exit\n";
			exit(1);
		}
		const char * bv2_vector = bv2.get_vector();
		for (unsigned long i = 0; i < boolean_vector_char_size; i++) {
			boolean_vector[i] = boolean_vector[i] &~ bv2_vector[i];
		}
	}
	
	void set_comment(const std::string & text)
	{
		comment = text;
	}
	
	const std::string & get_comment() const
	{
		return comment;
	}
};

#endif /* BOOLEAN_VECTOR_H_ */

