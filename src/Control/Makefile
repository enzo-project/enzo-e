
TOP     = ../..
PACKAGE = control
OBJS    = 
TESTS   = test_control
EXES    = cello

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

test_control: test_control.o 
	$(CXX) -o test_control test_control.o $(LIB)
#-----------------------------------------------------------------------
cello: cello.o 
	$(CXX) -o cello cello.o $(LIB)
#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/control.hpp      $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
