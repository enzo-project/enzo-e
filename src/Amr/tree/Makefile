
TOP = ../../..

EXES = test_tree

PG = 
INC_LOCAL  = -I$(TOP)/include
LIB_LOCAL  = $(LIB_HDF5)


OBJS = tree4.o node4.o tree16.o node16.o 


FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3 $(PG)
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3 $(PG)

FFLAGS   = $(FFLAGS_LOCAL)
CXXFLAGS = $(CXXFLAGS_LOCAL) $(INC_LOCAL)

all: $(EXES)
test_tree: $(OBJS) test_tree.o
	$(CXX) -o test_tree test_tree.o $(CXXFLAGS_LOCAL) $(OBJS) $(LIB)

obj: $(OBJS)

test: $(EXES)
#	./test_implosion # slow
#	./test_enzo      # slow

.PHONY: clean-local
clean-local:

include $(TOP)/src/MAKE/Make.include

include DEPEND
