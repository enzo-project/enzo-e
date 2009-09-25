
TOP     = ../..
PACKAGE = memory
OBJS    = memory.o
TESTS   = test_memory

include $(TOP)/src/MAKE/Make.include

#-----------------------------------------------------------------------

.PHONY: run-tests
run-tests: $(TESTS)
	-./test_memory

#-----------------------------------------------------------------------

test_memory: test_memory.o
	$(CXX) -o test_memory test_memory.o $(LIB)

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/memory.hpp $(INCLUDE)
