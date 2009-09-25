
TOP     = ../..
PACKAGE = parallel
OBJS    = parallel_parallel.o
TESTS   = 

include $(TOP)/src/MAKE/Make.include

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/parallel.hpp               $(INCLUDE)
	ln -sf `pwd`/parallel_parallel.hpp         $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
