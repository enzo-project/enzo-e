#======================================================================
#
# File:        Makefile
#
# Package:     amrSolve
#
# Description: Makefile for the amrSolve package
#
#----------------------------------------------------------------------
#
# Author:      James Bordner (jbordner@cosmos.ucsd.edu)
#
# History:     2002-12-03  created
#
#----------------------------------------------------------------------
#
# Copyright 2004 James Bordner
# Copyright 2004 Laboratory for Computational Astrophysics
# Copyright 2004 Regents of the University of California
#
#======================================================================

# User-configuration

ENZO_HOME = $(HOME)/enzo/svn/Enzo-devel/amr_mpi/src
ENZO_FLAGS = -D p4

OPTIMIZE = -O2

#----------------------------------------------------------------------

# Object files

# OBJ_LIB = grid.o

#----------------------------------------------------------------------

# Compiler flags

CXXFLAGS = $(ENZO_FLAGS) $(OPTIMIZE)

#----------------------------------------------------------------------

# Build targets

.PHONY: amr_solve library install

install: library

# library

library: $(OBJ_LIB)
	$(AR) -cr libamrSolve.a $(OBJ_LIB)
	mv libamrSolve.a ../lib

# main test program

amr_solve: $(OBJ_MAIN) $(OBJ_LIB)
	$(CXX) -o amr_solve $(OBJ_MAIN) $(OBJ_LIB)

# create enzo link

test-deriv:
	g++ -o test-discret -g test-discret.C

.PHONY: print

#=======================================================================
# print
#=======================================================================

CFILES = \
	grid.C  \
	hierarchy.C  \
	iterator.C   \
	main.C    \
	matrix.C  \
	print_hierarchy.C  \
	solver.C \
	vector.C

HFILES = \
	amrSolve.h \
	ENZO.h \
	grid.h \
	hierarchy.h \
	iterator.h \
	lGrid.h \
	main.h \
	matrix.h \
	MPI.h \
	solver.h \
	vector.h

print: $(patsubst %.C,print-C-%,$(CFILES)) $(patsubst %.h,print-h-%,$(HFILES))

print-h-%: 
	rm -f $*.hpp
	ln -s $*.h $*.hpp
	PRINT $*.hpp 1
	rm $*.hpp
	mv $*.hpp.ps $*.h.ps
	mv $*.h.ps ps

print-C-%: 
	PRINT $*.C 1
	mv $*.C.ps ps

#----------------------------------------------------------------------

# Other targets

clean:
	-rm -f *.o
	-rm -f *~
	-rm -f amr_solve
	-rm -f libamrSolve.a

