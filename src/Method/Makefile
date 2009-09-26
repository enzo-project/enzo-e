
TOP     = ../..
PACKAGE = method
OBJS    = method_method.o \
          method_ppm.o \
          method_ppm_solve_hydro_equations.o
TESTS   = 
DIRS    = ppm ppml

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

.PHONY: all
all: $(OBJS)
	$(MAKE) -C ppm
	$(MAKE) -C ppml

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/method.hpp               $(INCLUDE)
	ln -sf `pwd`/method_method.hpp         $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	-@(rm -f ppm/*.hdf5 ppm/*.png ppm/*.mpg)

