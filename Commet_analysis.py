# Contributors :
#   Pierre PETERLONGO, pierre.peterlongo@inria.fr [12/06/13]
#   Nicolas MAILLET, nicolas.maillet@inria.fr     [12/06/13]
#   Guillaume Collet, guillaume@gcollet.fr        [27/05/14]
#
# This software is a computer program whose purpose is to find all the
# similar reads between sets of NGS reads. It also provide a similarity
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
from Commet import getReadBVFiles
from Commet import fillDefaultBVReadSetMatrix

def main():
    parser = argparse.ArgumentParser(description='Computes the matrices from .bv results')
    parser.add_argument("input_file", type=str,
                        help="input file of files (a line=a set composed by: \"set_name:read_file;read_file;read_file...\")" )
                    

    
    parser.add_argument("-b", "--binaries_directory", type=str, dest='binary_directory', metavar='',
                        help="binary directory  [default: \"./bin\"]", default="./bin" )
                    
    parser.add_argument("-o", "--output_directory", type=str, dest='directory', metavar='',
                        help="directory in which vector results will be output [default: \"output_commet\"]", default="output_commet/" )
  
    

    args = parser.parse_args()
 
    # The input file of files
    input_file=str(args.input_file)
    output_directory=str(args.directory)
    if output_directory[-1]!='/': output_directory+="/"
    bin_dir=str(args.binary_directory)
    if bin_dir[-1]!='/': bin_dir+="/"


    # Stores the input reads in a matrix
    readSetMatrix = getReadFiles(input_file)
    bvreadSetMatrix = getReadBVFiles(input_file)
    if bvreadSetMatrix == None:
        bvreadSetMatrix = fillDefaultBVReadSetMatrix(readSetMatrix, output_directory)
    
    readSetNames = getReadSetsNames(input_file)
    
    output_matrices (readSetMatrix, bvreadSetMatrix, readSetNames, output_directory, bin_dir)
    

if __name__ == "__main__":
    main()
        
 











