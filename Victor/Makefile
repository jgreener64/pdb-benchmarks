# Makefile for Victor benchmarks
# VICTOR_ROOT needs to be set as the Victor root directory

BINPATH = $(VICTOR_ROOT)/bin

LIBS = -lBiopool -ltools
LIB_PATH = -L. -L$(VICTOR_ROOT)/lib/
INC_PATH += -I. -I$(VICTOR_ROOT)/tools/ -I$(VICTOR_ROOT)/Biopool/Sources/

CC=g++
CFLAGS=-I. -ansi -pedantic -DNEXCEPTIONS -DLINUX -c -O3 -ffast-math -DNDEBUG -ftemplate-depth-36 -Wno-reorder  -Wno-uninitialized -Wno-write-strings -Wno-narrowing

install:
	$(CC)   $(CFLAGS)   $(INC_PATH) -c parse_pdb.cc -o parse_pdb.o
	$(CC)   parse_pdb.o -o parse_pdb $(LIB_PATH) $(LIBS)

clean:
	rm parse_pdb.o parse_pdb
