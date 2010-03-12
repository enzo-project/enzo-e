TOP     = ../..
PACKAGE = array
OBJS    = array.o
TESTS   = test_array

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

test_array: test_array.o 
	$(CXX) -o test_array test_array.o $(LIB)
#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/array.hpp $(INCLUDE)
	ln -sf `pwd`/array_array.hpp $(INCLUDE)
	ln -sf `pwd`/array_block.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
