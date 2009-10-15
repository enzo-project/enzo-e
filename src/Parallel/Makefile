
TOP     = ../..
PACKAGE = parallel
OBJS    = parallel.o mpi.o
TESTS   = 

include $(TOP)/src/MAKE/Make.include

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/parallel.hpp               $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
