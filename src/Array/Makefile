TOP     = ../..
PACKAGE = array
OBJS    = array_serial.o
TESTS   = test_array

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
	-./test_array
#-----------------------------------------------------------------------
test_array: test_array.o 
	$(CXX) -o test_array test_array.o $(LIB)
#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/array.hpp         $(INCLUDE)
	ln -sf `pwd`/array_serial.hpp  $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
