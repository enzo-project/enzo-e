# (X) install:         Install lib, include, and bin files
# (X) install-install: Install lib, include, and bin files
# 

TOP = ../..
SRC = ..
LIB = -L../../lib -lcello

TEST_EXES = \
	test_array \
	test_disk_hdf5 \
	test_error \
	test_parameters \
	test_performance_timer

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: test run-tests

#-----------------------------------------------------------------------

.PHONY: test
test: $(TEST_EXES)

#-----------------------------------------------------------------------

.PHONY: run-tests
run-tests: $(TEST_EXES)
	-./test_array
	-./test_disk_hdf5
	-./test_error
	-./test_parameters
	-./test_performance_timer
#-----------------------------------------------------------------------

test_array: test_array.o 
	$(CXX) -o test_array test_array.o $(LIB) $(INC)
test_disk_hdf5: test_disk_hdf5.o
	$(CXX) -o test_disk_hdf5 test_disk_hdf5.o $(LIB) $(LIB_HDF5) $(INC)
test_error: test_error.o
	$(CXX) -o test_error test_error.o $(LIB) $(INC)
test_parameters: test_parameters.o
	$(CXX) -o test_parameters test_parameters.o $(LIB) $(INC)
test_performance_timer: test_performance_timer.o 
	$(CXX) -o test_performance_timer test_performance_timer.o $(LIB) $(INC)


test_test:

#-----------------------------------------------------------------------

libtest.a:

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ $(TEST_EXES) *.unit *.hdf5

#-----------------------------------------------------------------------

.PHONY: install
install: $(TEST_EXES)
	cp $(TEST_EXES) $(TOP)/bin


#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/test.hpp         $(INCLUDE)
	ln -sf `pwd`/test_unit_.hpp   $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include
