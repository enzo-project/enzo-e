
TOP     = ../..
PACKAGE = schedule
OBJS    = schedule.o
TESTS   = test_schedule
EXES    = cello

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

test_schedule: test_schedule.o 
	$(CXX) -o test_schedule test_schedule.o $(LIB)
#-----------------------------------------------------------------------
cello: cello.o 
	$(CXX) -o cello cello.o $(LIB)
#-----------------------------------------------------------------------

.PHONY: install-include
install-include:

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
