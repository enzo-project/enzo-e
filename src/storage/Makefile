#
# (x) test:         Run unit tests
# (x) test_array:   Compile the unit tests
# (x) libarray.a:
# (x) clean:        Clean directory
# (^) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

TEST_OBJS = test_hdf5.o 
OBJ       = hdf5.o

#=======================================================================
INC_LOCAL  = -I$(SRC)/data -I$(SRC)/error
OBJ_EXTERN = $(SRC)/data/array.o

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================

#=======================================================================
# TARGETS
#=======================================================================

#-----------------------------------------------------------------------

.PHONY: test
test: test_hdf5
	./test_hdf5

#-----------------------------------------------------------------------

test_hdf5: $(TEST_OBJS) $(OBJS)
	$(CXX) -o test_hdf5 $(TEST_OBJS) $(OBJS) $(INC) $(LIB_HDF5)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-rm *.o test_*hdf5 *~ out.make.dep Makefile.bak lib*.a *.unit *.hdf5

#-----------------------------------------------------------------------

.PHONY: install
install: install-include test_hdf5 $(OBJS)
	cp ./test_hdf5 $(TOP)/bin
	ar rcu $(TOP)/lib/libcello.a $(OBJS)

#-----------------------------------------------------------------------

.PHONY: install-include

install-include:
	ln -sf `pwd`/storage.hpp  $(TOP)/include
	ln -sf `pwd`/hdf5.hpp     $(TOP)/include/storage_hdf5_.hpp

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include




