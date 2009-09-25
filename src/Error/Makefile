#
# ( ) liberror.a:
# ( ) clean:        Clean directory
# ( ) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..
LIB = -L$(TOP)/lib -lcello

PACKAGE  = error
OBJS     = 
TESTS    = test_error

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

.PHONY: clean
clean:
	-@rm *.o *~ lib*.a *.unit

#-----------------------------------------------------------------------

.PHONY: install

install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/error.hpp            $(INCLUDE)
	ln -sf `pwd`/error_exception.hpp  $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include




