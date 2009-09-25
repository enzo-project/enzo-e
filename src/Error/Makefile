
TOP     = ../..
PACKAGE = error
OBJS    = 
TESTS   = test_error

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: $(OBJS) $(TESTS)

#-----------------------------------------------------------------------

.PHONY: test
test: $(TEST_EXES) run-tests

#-----------------------------------------------------------------------

.PHONY: run-tests
run-tests: $(TESTS)
	-./test_error

#-----------------------------------------------------------------------
test_error: test_error.o 
	$(CXX) -o test_error test_error.o $(LIB) $(LIB_HDF5)

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/error.hpp            $(INCLUDE)
	ln -sf `pwd`/error_exception.hpp  $(INCLUDE)
