
TOP     = ../..
PACKAGE = performance
OBJS    = 
TESTS   = test_performance
EXES    = time_matvec

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

test_performance: test_performance.o
	$(CXX) -o test_performance test_performance.o $(LIB)

#-----------------------------------------------------------------------

# force recompile of time_matvec in case OPT changes
.PHONY: time_matvec
time_matvec: time_matvec.cpp
	@echo "Compile time_matvec with, e.g., OPT=-O3"
	$(CXX) -o time_matvec.o time_matvec.cpp $(OPT) $(INC) -c
	$(CXX) -o $@ time_matvec.o  $(OPT) $(LIB)

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/performance.hpp       $(INCLUDE)
	ln -sf `pwd`/performance_timer.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:

