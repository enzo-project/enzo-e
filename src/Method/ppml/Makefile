#
# test         Run unit tests
# clean        Clean directory
# test_ppml    Compile the unit tests
# dep          Generate dependencies
# 

TOP = ../../..

PACKAGE = method

TESTS = test_ppml

INC_LOCAL  = 
LIB_LOCAL  = -L$(TOP)/lib -lcello $(LIB_HDF5) $(LIB_FORTRAN)
CPPFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3

INCS = \

OBJS = \
	cello_hydro.o \
	../CosmologyComputeExpansionFactor.o \
	../CosmologyComputeExpansionTimestep.o \
	../data_dump.o \
        Grid_ComputeTimeStep.o \
	../Grid_SetExternalBoundaryValues.o \
	Grid_SolveMHDEquations.o \
	../initialize_hydro.o \
        initialize_ppml.o \
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
	../SetExternalBoundary.o

include $(TOP)/src/MAKE/Make.include

.PHONY: clean-local
clean-local:
	-@rm -f *.png *.hdf5


test_ppml: $(INCS) $(OBJS) test_ppml.o
	$(CXX) -o test_ppml test_ppml.o $(OBJS) $(LIB)


.PHONY: includes
includes:
	ln -s ../ppm/cello_hydro.h .
	ln -s ../ppm/defines.h .


include DEPEND
