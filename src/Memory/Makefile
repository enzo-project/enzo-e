#
# (x) clean:        Clean directory
# (^) dep:          Update dependencies
# (X) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..
LIB = -L$(TOP)/lib -lcello

PACKAGE = memory
OBJS    = memory.o
TESTS   = test_memory

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
	-./test_memory

#-----------------------------------------------------------------------

test_memory: test_memory.o
	$(CXX) -o test_memory test_memory.o $(LIB)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ out.make.dep Makefile.bak $(TESTS)

#-----------------------------------------------------------------------

.PHONY: install
install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)   
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE = $(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/memory.hpp $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include

# DO NOT DELETE
