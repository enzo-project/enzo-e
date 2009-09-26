
TOP     = ../..
PACKAGE = memory
OBJS    = memory.o
TESTS   = test_memory

include $(TOP)/src/MAKE/Make.include

#-----------------------------------------------------------------------

test_memory: test_memory.o
	$(CXX) -o test_memory test_memory.o $(LIB)

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/memory.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
