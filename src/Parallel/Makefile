
TOP     = ../..
PACKAGE = parallel
OBJS    = parallel.o parallel_mpi.o
TESTS   = 
DIRS    = mpi-tests omp-tests

include $(TOP)/src/MAKE/Make.include

#=======================================================================

# override all target to compile $(DIRS)

 .PHONY: all
all: $(OBJS) $(TESTS) $(EXES)
	$(MAKE) -C mpi-tests
	$(MAKE) -C omp-tests

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/parallel.hpp $(INCLUDE)
	ln -sf `pwd`/parallel_parallel.hpp $(INCLUDE)
	ln -sf `pwd`/parallel_mpi.hpp $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	$(MAKE) -C mpi-tests clean
	$(MAKE) -C omp-tests clean
