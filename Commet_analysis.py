# Contributors :
#   Pierre PETERLONGO, pierre.peterlongo@inria.fr [12/06/13]
#   Nicolas MAILLET, nicolas.maillet@inria.fr     [12/06/13]
#   Guillaume Collet, guillaume@gcollet.fr        [27/05/14]
#
# This software is a computer program whose purpose is to find all the
# similar reads between two set of NGS reads. It also provide a similarity
# score between the two samples.
#
# Copyright (C) 2014  INRIA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import errno
import random
import string
import argparse
import subprocess
from Commet import getReadFiles
from Commet import getReadSetsNames
from Commet import output_matrices



def main():
    parser = argparse.ArgumentParser(description='Computes the matrices from .bv results')
    parser.add_argument("input_file", type=str,
                        help="input file of files (a line=a set composed by: \"set_name:read_file;read_file;read_file...\")" )
                    
    parser.add_argument("-o", "--output_directory", type=str, dest='directory', metavar='',
                        help="directory in which vector results will be output [default: \"output_commet\"]", default="output_commet/" )
  
    parser.add_argument("prefix_matrix", type=str,
                        help="prefix of files in which matrices are output (in the output directory)[default=stdout]", default=None)
    

    args = parser.parse_args()
 
    # The input file of files
    input_file=str(args.input_file)
    output_directory=str(args.directory)+"/"
    output_matrix_prefix=args.prefix_matrix
    print "input file="+input_file, 
    

    #ouput directory
    print " output directory="+output_directory

    print "output matrices in:"
    print "\t"+output_directory+output_matrix_prefix+"_plain.csv"
    print "\t"+output_directory+output_matrix_prefix+"_percentage.csv"
    print "\t"+output_directory+output_matrix_prefix+"_normalized.csv"


    # Stores the input reads in a matrix
    readSetMatrix = getReadFiles(input_file)
    readSetNames = getReadSetsNames(input_file)
    output_matrices (readSetMatrix, readSetNames, output_directory, output_matrix_prefix)
    
    print "All Comme_analysis work is done, output matrices are in:"
    print "\t"+output_directory+output_matrix_prefix+"_plain.csv and .pdf"
    print "\t"+output_directory+output_matrix_prefix+"_percentage.csv and .pdf"
    print "\t"+output_directory+output_matrix_prefix+"_normalized.csv and .pdf"

if __name__ == "__main__":
    main()
        
 











