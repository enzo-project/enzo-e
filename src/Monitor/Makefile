
TOP     = ../..
PACKAGE = monitor
OBJS    = 
TESTS   = 

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

# test_monitor: test_monitor.o
# 	$(CXX) -o test_monitor test_monitor.o $(LIB)

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/monitor.hpp       $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:

