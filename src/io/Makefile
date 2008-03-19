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
OBJ        = hdf5.o

#=======================================================================
INC_LOCAL  = -I$(SRC)/data
OBJ_EXTERN = $(SRC)/data/array.o

OBJS = $(OBJ) $(OBJ_EXTERN)
#=======================================================================

#=======================================================================
# TARGETS
#=======================================================================

#-----------------------------------------------------------------------

.PHONY: test
test: test_hdf5
	test_hdf5

#-----------------------------------------------------------------------

test_hdf5: $(TEST_OBJS) libio.a
	$(CPP) -o test_hdf5 $(TEST_OBJS) $(INC) $(INC_IO) -lhdf5 -L. -lio

#-----------------------------------------------------------------------

libio.a: $(OBJS)
	ar rv $*.a $(OBJS)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-rm *.o test_*hdf5 *~ out.make.dep Makefile.bak lib*.a *.unit

#-----------------------------------------------------------------------

.PHONY: install
install: test_hdf5 libio.a
	cp ./test_hdf5 $(TOP)/bin
	ar rv $(TOP)/lib/libcello.a libio.a

#-----------------------------------------------------------------------

.PHONY: install-include

install-include:
	ln -sf `pwd`/io.hpp    $(TOP)/include
	ln -sf `pwd`/hdf5.hpp  $(TOP)/include/io_hdf5_.hpp
	ln -sf `pwd`/error.hpp $(TOP)/include/io_error_.hpp

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/Make.include




