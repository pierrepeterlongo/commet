/*
 * Contributors :
 *   Pierre PETERLONGO, pierre.peterlongo@inria.fr [09/07/18]
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

#include "boolean_vector.h"
#include "file_manager.h"
#include <iostream>

std::string version = "2.1";

// -----------------------------------------------------------------------
//                             PRINT USAGE
// -----------------------------------------------------------------------
void print_usage ()
{
	std::cout << "\ngenerate_random_bv, version " << version << "\n";
	std::cout << "Usage : ./generate_random_bv <read_set> <percentage_kept_reads> <output_bv_name>\n";
	std::cout << "Mandatory:\n";
	std::cout << "\t<read_set>\t: read file \n";
	std::cout << "\t<percentage_kept_reads>\t: percentage of reads to be kept from read_set \n";
	std::cout << "\t<output_bv_name>\t: name of the resulting boolean vector file name \n";
}


int main (int argc, char ** argv)
{
	
	////////////////////////////////////////////////////////////
	// Read command line arguments through getopt
	//
	if (argc < 4) {
		std::cerr << "A a read file, an integer and an ouput bv name must be provided. See usage.\n";
		print_usage();
		return(1);
	}
	std::string read_set = argv[1];
    const int percentage_kept_reads = atoi(argv[2]);
    if (percentage_kept_reads<0 or percentage_kept_reads>100){
        std::cerr <<" the Percentage of read to be kept must be in [0,100]. See usage.\n";
		print_usage();
		return(1);
	}
    
	std::string output_bv_name = argv[3];
    
    FileManager fm =  FileManager();
    fm.addFile(read_set);
    std::vector<BooleanVector> bvs = fm.get_bvs ();
    
    BooleanVector bv = bvs[0]; 
    bv.random_vector(percentage_kept_reads);
    std::stringstream comment;
    comment << percentage_kept_reads << " % random reads kept";
	bv.set_comment(comment.str());
	bv.print(output_bv_name);

	return 0;
}
