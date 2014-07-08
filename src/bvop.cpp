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

#include "boolean_vector.h"
#include "file_manager.h"
#include <iostream>

std::string version = "2.1";

// -----------------------------------------------------------------------
//                             PRINT USAGE
// -----------------------------------------------------------------------
void print_usage ()
{
	std::cout << "\nbvop, version " << version << "\n";
	std::cout << "Usage : ./bvop <file1.bv> [options]\n";
	std::cout << "Mandatory:\n";
	std::cout << "\t<file1.bv>\t: file containing a boolean vector\n";
	std::cout << "Options:\n";
	std::cout << "\t -n             : performs NOT on file1.bv\n";
	std::cout << "\t -a <file2.bv>  : performs file1.bv AND file2.bv\n";
	std::cout << "\t -o <file2.bv>  : performs file1.bv OR file2.bv\n";
	std::cout << "\t -d <file2.bv>  : performs file1.bv AND (NOT file2.bv)\n";
	std::cout << "\t -p <output.bv> : print result in file output.bv [Default=stdout]\n";
	std::cout << "\t -i             : print information about file1.bv\n";
	std::cout << "\t -h             : Prints this message and exit\n";
	std::cout << "\t -v             : Prints the version number and exit\n";
}


int main (int argc, char ** argv)
{
	
	////////////////////////////////////////////////////////////
	// Read command line arguments through getopt
	//
	if (argc < 2) {
		std::cerr << "A boolean vector file must be provided, see usage\n";
		print_usage();
		return(1);
	}
	std::string file_name1;
	std::string file_name2;
	std::string output_file_name;
	bool print = false;
	bool print_info = false;
	bool do_nothing = false;
	char bvop = 'u';
	int i = 1;
	while (i < argc)
	{
		if (argv[i][0] == '-') {
			int flag = argv[i][1];
			switch (flag)
			{
				case 'a':
					i++;
					file_name2 = argv[i];
					bvop = 'a';
					break;
				case 'o':
					i++;
					file_name2 = argv[i];
					bvop = 'o';
					break;
				case 'd':
					i++;
					file_name2 = argv[i];
					bvop = 'd';
					break;
				case 'n':
					bvop = 'n';
					break;
				case 'p':
					i++;
					output_file_name = argv[i];
					print = true;
					break;
				case 'i':
					print_info = true;
					break;
				case 'v':
					// Print the version
					std::cout << "compare_reads version " << version << "\n";
					return (0);
				case 'h':
					// Print the help
				default:
					print_usage ();
					return (0);
			}
		} else {
			if (file_name1.empty()) {
				file_name1 = argv[i];
			} else {
				std::cerr << "One input file is mandatory\n";
				print_usage ();
				return (0);
			}
		}
		i++;
	}
	
	////////////////////////////////////////////////////////////
	// Read the boolean vectors from files
	//
	BooleanVector bv1;
	bv1.read(file_name1);
	std::string comment;
	if (bvop == 'a') {
		BooleanVector bv2;
		bv2.read(file_name2);
		bv1.full_and(bv2);
		comment = file_name1 + " AND " + file_name2 + "\n";
	} else if (bvop == 'o') {
		BooleanVector bv2;
		bv2.read(file_name2);
		bv1.full_or(bv2);
		comment = file_name1 + " OR " + file_name2 + "\n";
	} else if (bvop == 'd') {
		BooleanVector bv2;
		bv2.read(file_name2);
		bv1.full_and_not(bv2);
		comment = file_name1 + " AND (NOT " + file_name2 + ")\n";
	} else if (bvop == 'n') {
		bv1.full_not();
		comment = "NOT " + file_name1 + "\n";
	} else {
		do_nothing = true;
	}
	
	if (print_info) {

		std::cout << bv1.get_comment();
		std::cout << "\nReads:\n";
		std::cout << "  " << bv1.nb_one() << " / " << bv1.size() << " reads selected\n";
	}
	
	if (do_nothing) {
		return 0;
	}
	
	bv1.set_comment(comment);
	if (print) {
		bv1.print(output_file_name);
	} else {
		bv1.print();
	}
	
	
	return 0;
}
