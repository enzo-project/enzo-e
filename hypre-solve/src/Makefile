SHELL = /bin/bash

# Define MPCPP
include Make.machine

FLAGS = -g

# JBPERF

JBLIB = -ljbperf
JBDEF = -DUSE_JBPERF

# Assemble defines, includes, and libs
DEF = $(JBDEF)
INC = -I$(HOME)/include 
LIBS = -L$(HOME)/lib -lHYPRE_sstruct_ls -lHYPRE_sstruct_mv -lHYPRE_struct_ls -lHYPRE_struct_mv -lHYPRE_parcsr_ls -lHYPRE_parcsr_mv -lHYPRE_IJ_mv -lHYPRE_LSI -lHYPRE_MatrixMatrix -lHYPRE_DistributedMatrix -lHYPRE_utilities $(JBLIB)

# Set object files
OBJS = \
       domain.o \
       faces.o \
       grid.o \
       hierarchy.o \
       hypre.o \
       hypre-solve.o \
       level.o \
       mpi.o \
       parameters.o \
       point.o \
       problem.o \
       sphere.o

hypre-solve: dep $(OBJS)
	$(MPCPP) -o hypre-solve $(OBJS) $(LIBS)

BIN = ../bin
.PHONY: install
install: hypre-solve
	(if [ ! -e $(BIN) ]; then mkdir $(BIN); fi)
	cp hypre-solve $(BIN)


.cpp.o:
	$(MPCPP) -c -o $*.o $*.cpp $(FLAGS) $(INC) $(DEF)

.PHONY: dep
dep:
	-@(makedepend -fDEPEND *.cpp *.hpp) >& /dev/null

include DEPEND

.PHONY: clean
clean:
	rm -f *.o hypre-solve *~ \#*

.PHONY: distclean
distclean: clean
	rm -rf html latex

.PHONY: doc
doc:
	doxygen


