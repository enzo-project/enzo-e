#
# (x) test:         Run unit tests
# (x) test_array:   Compile the unit tests
# (x) libdata.a:
# (x) clean:        Clean directory
# (^) dep:          Update dependencies
# (X) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..

TEST_OBJS = test_array.o 
OBJS = array.o

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: test
test: test_array
	test_array

#-----------------------------------------------------------------------

test_array: $(TEST_OBJS) libdata.a
	$(CPP) -o test_array $(TEST_OBJS) $(INC) -L. -ldata

#-----------------------------------------------------------------------

libdata.a: $(OBJS)
	ar rv $*.a $(OBJS)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-rm *.o test_array *~ out.make.dep Makefile.bak lib*.a *.unit

#-----------------------------------------------------------------------

.PHONY: install
install: test_array libdata.a install-include
	cp ./test_array $(TOP)/bin
	ar rv $(TOP)/lib/libcello.a libdata.a   

#-----------------------------------------------------------------------

.PHONY: install-include

install-include:
	ln -sf `pwd`/data.hpp   $(TOP)/include
	ln -sf `pwd`/array.hpp  $(TOP)/include/data_array_.hpp
	ln -sf `pwd`/scalar.hpp $(TOP)/include/data_scalar_.hpp

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/Make.include

