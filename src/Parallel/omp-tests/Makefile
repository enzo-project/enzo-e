TOP     = ../../..
PACKAGE = parallel
OBJS    = omp-test.o
TESTS   = omp-test

CXXFLAGS_LOCAL = -fopenmp
LIB_LOCAL      = -fopenmp

include $(TOP)/src/MAKE/Make.include

#-----------------------------------------------------------------------

omp-test: omp-test.o
	$(CXX) -o omp-test omp-test.o $(LIB) $(LIB_HDF5)

.PHONY: install-include
install-include:

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
