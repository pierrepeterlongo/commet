# Contributors :
#   Pierre PETERLONGO, pierre.peterlongo@inria.fr [12/06/13]
#   Nicolas MAILLET, nicolas.maillet@inria.fr     [12/06/13]
#   Guillaume Collet, guillaume@gcollet.fr        [22/01/14]
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
#

CFLAGS=-O3 -Wall -Iinclude/
LDFLAGS=-lm -lz
CC=g++
HDRS= $(wildcard include/*.h)

ifeq ($(prof),1)
    CFLAGS= -O3 -pg -g
endif


all: bin/index_and_search bin/filter_reads bin/extract_reads bin/bvop bin/compare_reads

bin/index_and_search: src/index_and_search.cpp $(HDRS)
	@ if [ ! -d bin ]; then mkdir bin; fi
	$(CC) -o bin/index_and_search src/index_and_search.cpp $(LDFLAGS) $(CFLAGS)

bin/filter_reads: src/filter_reads.cpp $(HDRS)
	@ if [ ! -d bin ]; then mkdir bin; fi
	$(CC) -o bin/filter_reads src/filter_reads.cpp $(LDFLAGS) $(CFLAGS)

bin/compare_reads: src/compare_reads.cpp $(HDRS)
	@ if [ ! -d bin ]; then mkdir bin; fi
	$(CC) -o bin/compare_reads src/compare_reads.cpp $(LDFLAGS) $(CFLAGS)

bin/extract_reads: src/extract_reads.cpp $(HDRS)
	@ if [ ! -d bin ]; then mkdir bin; fi
	$(CC) -o bin/extract_reads src/extract_reads.cpp $(LDFLAGS) $(CFLAGS)

bin/bvop: src/bvop.cpp $(HDRS)
	@ if [ ! -d bin ]; then mkdir bin; fi
	$(CC) -o bin/bvop src/bvop.cpp $(LDFLAGS) $(CFLAGS)

install:
	cp bin/filter_reads /usr/local/bin/
	cp bin/extract_reads /usr/local/bin/
	cp bin/bvop /usr/local/bin/
	cp bin/index_and_search /usr/local/bin/
clean:
	@ rm bin/*
