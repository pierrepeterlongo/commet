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

#ifndef __READ_FILE_H__
#define __READ_FILE_H__

#include "boolean_vector.h"

//
// Interface to files handling reads
//
class ReadFile
{
protected:
	std::string fname;
	std::string current_read_data;
	std::string current_read_seq;
	unsigned long _nb_valid_reads;
	unsigned long _cnt_valid_reads;
	unsigned long current_read_pos;
	unsigned long nb_reads;
	BooleanVector bv;
	bool first_read;
public:
	virtual ~ReadFile () {};
	virtual std::string & get_next_read () = 0;
	virtual const std::string & get_data() const = 0;
	virtual void flush_next_read () = 0;
	virtual void tag_current_read () = 0;
	virtual void untag_current_read () = 0;
	virtual void rewind () = 0;
	virtual void save_bv (const std::string & file_name) = 0;
	virtual void save_bv () = 0;
	virtual void save (const std::string & directory, const std::string & suffix) = 0;
	virtual const unsigned long & get_nb_reads() const {return nb_reads;}
	virtual const unsigned long & get_read_pos() const {return current_read_pos;}
	virtual void tag (const unsigned long & pos) = 0;
	virtual void untag (const unsigned long & pos) = 0;
	virtual void set_bv_comment (const std::string & str) = 0;
	virtual const BooleanVector & get_bv () const {return bv;}
	virtual const BooleanVector * get_bv_ptr () const {return &bv;}
	virtual unsigned long nb_valid_reads () {return _nb_valid_reads;}
	virtual const std::string & get_fname () const {return fname;}
	virtual void apply_bv(const BooleanVector & new_bv)
	{
		bv = new_bv;
		_nb_valid_reads = bv.nb_one();
	};
	
	////////////////////////////////////////////////////////////
	// Unset the bits of the reads from beg to beg+nb
	//
	void untag_last_reads ()
	{
		for (unsigned long i = current_read_pos; i < nb_reads; i++) {
			bv.unset(i);
		}
	}

};

#endif
