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
FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3

INCS = \
	cello_hydro.h \
	data_dump.cpp \
	SetExternalBoundary.cpp \
	initialize_hydro.cpp \
	CosmologyComputeExpansionTimestep.cpp \
	Grid_SetExternalBoundaryValues.cpp \
	CosmologyComputeExpansionFactor.cpp \

OBJS = \
	cello_hydro.o \
        initialize_hydro.o \
        initialize_ppml.o \
	CosmologyComputeExpansionFactor.o \
	CosmologyComputeExpansionTimestep.o \
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
	data_dump.o


# all: $(OBJS)
#	PPML_QDD6.o \

include $(TOP)/src/MAKE/Make.include

.PHONY: clean-local
clean-local:
	-@rm -f *.png *.hdf5


test_ppml: $(INCS) $(OBJS) test_ppml.o
	$(CXX) -o test_ppml test_ppml.o $(OBJS) $(LIB)

cello_hydro.h:
	ln -s ../ppm/cello_hydro.h .

data_dump.cpp:
	ln -s ../ppm/data_dump.cpp .

SetExternalBoundary.cpp:
	ln -s ../ppm/SetExternalBoundary.cpp .

initialize_hydro.cpp:
	ln -s ../ppm/initialize_hydro.cpp .

CosmologyComputeExpansionTimestep.cpp:
	ln -s ../ppm/CosmologyComputeExpansionTimestep.cpp .

Grid_SetExternalBoundaryValues.cpp:
	ln -s ../ppm/Grid_SetExternalBoundaryValues.cpp .

CosmologyComputeExpansionFactor.cpp:
	ln -s ../ppm/CosmologyComputeExpansionFactor.cpp .


include DEPEND
