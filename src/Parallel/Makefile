
TOP     = ../..
PACKAGE = parallel
OBJS    = parallel.o mpi.o
TESTS   = 
DIRS    = mpi-tests

include $(TOP)/src/MAKE/Make.include

#=======================================================================

# override all target to compile $(DIRS)

 .PHONY: all
all: $(OBJS) $(TESTS) $(EXES)
	$(MAKE) -C mpi-tests

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/parallel.hpp               $(INCLUDE)
	ln -sf `pwd`/mpi.hpp                    $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	$(MAKE) -C mpi-tests clean

