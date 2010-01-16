
TOP     = ../..
PACKAGE = memory
OBJS    = memory.o
TESTS   = test_memory

include $(TOP)/src/MAKE/Make.include

#-----------------------------------------------------------------------

test_memory: test_memory.o $(OBJS)
	$(CXX) -o test_memory test_memory.o $(OBJS)

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/memory.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
