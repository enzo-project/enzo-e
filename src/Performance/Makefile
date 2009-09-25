#
# ( ) libperformance.a:
# ( ) clean:        Clean directory
# ( ) dep:          Update dependencies
# ( ) install:      Install lib, include, and bin files
# 

TOP = ../..
SRC = ..
LIB = -L$(TOP)/lib -lcello

PACKAGE   = performance
OBJS      = 
TESTS     = test_performance time_matvec

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
	-./test_performance
	-./time_matvec

#-----------------------------------------------------------------------

test_performance: test_performance.o
	$(CXX) -o test_performance test_performance.o $(LIB)

#-----------------------------------------------------------------------

# force recompile of time_matvec in case OPT changes
.PHONY: time_matvec
time_matvec: time_matvec.cpp
	@echo "Compile time_matvec with, e.g., OPT=-O3"
	$(CXX) -o time_matvec.o time_matvec.cpp $(OPT) $(INC) -c
	$(CXX) -o $@ time_matvec.o  $(OPT) $(INC)

#-----------------------------------------------------------------------

.PHONY: clean
clean:
	-@rm *.o *~ $(TESTS) out.make.dep Makefile.bak lib*.a *.unit *.png *.eps

#-----------------------------------------------------------------------

.PHONY: install
install: install-include $(OBJS)
	ar rcu $(TOP)/lib/libcello.a      $(OBJS)
	ar rcu $(TOP)/lib/lib$(PACKAGE).a $(OBJS)   

#-----------------------------------------------------------------------

INCLUDE =$(TOP)/include

.PHONY: install-include
install-include:
	ln -sf `pwd`/performance.hpp       $(INCLUDE)
	ln -sf `pwd`/performance_timer.hpp $(INCLUDE)

#=======================================================================
# INCLUDES 
#=======================================================================

include DEPEND
include $(SRC)/MAKE/Make.include




