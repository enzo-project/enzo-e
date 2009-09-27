#
# test         Run unit tests
# clean        Clean directory
# test_ppml    Compile the unit tests
# dep          Generate dependencies
# 

TOP = ../../..
PACKAGE = method
OBJS = \
	cello_hydro.o \
        initialize_hydro.o \
        initialize_implosion.o \
	Grid_SolveMHDEquations.o \
	SetExternalBoundary.o \
        Grid_SetExternalBoundaryValues.o \
        Grid_ComputeTimeStep.o \
	PPML_Conservative.o \
	PPML_HLLD.o \
	PPML_MAIN.o \
	PPML_monot.o \
	PPML_potokx.o \
	PPML_potoky.o \
	PPML_potokz.o \
	PPML_Primitive.o \
	PPML_sub.o \
	PPML_TimeStep.o \
	calc_dt.o 
TESTS = test_ppml
all: $(OBJS)
#	PPML_QDD6.o \

include $(TOP)/src/MAKE/Make.include

FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3

LIB      = $(LIB_LOCAL) -lgfortran -lm -lhdf5

.PHONY: clean-local
clean-local:
	-@rm -f *.png *.hdf5


test_ppml: $(OBJS) test_ppml.o
	$(CXX) -o test_ppml test_ppml.o $(OBJS) $(LIB)


include DEPEND
