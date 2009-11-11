
TOP     = ../..
PACKAGE = method
OBJS    = method.o \
          method_ppm.o \
          method_ppm_solve_hydro_equations.o
TESTS   = 
DIRS    = ppm ppml

include $(TOP)/src/MAKE/Make.include

#=======================================================================
# TARGETS
#=======================================================================

# override all target to compile ppm and ppml
 .PHONY: all
all: $(OBJS) $(TESTS) $(EXES)
	$(MAKE) -C ppm
	$(MAKE) -C ppml

#-----------------------------------------------------------------------

.PHONY: install-include
install-include:
	ln -sf `pwd`/defines.h                 $(INCLUDE)
	ln -sf `pwd`/method.hpp                $(INCLUDE)
	ln -sf `pwd`/method_ppm.hpp                $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	-@$(MAKE) -C ppm clean
	-@$(MAKE) -C ppml clean
	-@(rm -f ppm/*.hdf5 ppm/*.png ppm/*.mpg)

