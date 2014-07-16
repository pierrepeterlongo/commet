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

#ifndef ALPHABET_H_
#define ALPHABET_H_

/*
 * Alphabet is a singleton class (instanciate once and only once)
 *
 * It is used to know if a char is an A, C, G, or T
 * See usage in index_reads.h and search_reads.h
 */
class Alphabet
{
public:
	static Alphabet * getInstance()
	{
		static Alphabet theOnlyInstance;
		return &theOnlyInstance;
	}
	bool is_in (int i) {return (bool) alpha[i];};
private:
	Alphabet() {
		for (int i = 0; i < 255; i++) {
			alpha[i] = 0;
		}
		alpha[(int)'A'] = 1;
		alpha[(int)'a'] = 1;
		alpha[(int)'C'] = 1;
		alpha[(int)'c'] = 1;
		alpha[(int)'G'] = 1;
		alpha[(int)'g'] = 1;
		alpha[(int)'T'] = 1;
		alpha[(int)'t'] = 1;
	};
	~Alphabet() {};
	
	char alpha [255];
};

#endif
