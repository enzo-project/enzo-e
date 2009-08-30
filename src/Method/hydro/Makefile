#
# test         Run unit tests
# clean        Clean directory
# test_array   Compile the unit tests
# dep          Generate dependencies
# 

TOP = ../../..

INC_LOCAL  = -I$(TOP)/include

OBJS = \
	CosmologyComputeExpansionFactor.o \
	CosmologyComputeExpansionTimestep.o \
	CosmologyGetUnits.o \
        cello_hydro.o \
        c_message.o \
	FindField.o \
	Grid_ComputeGammaField.o \
	Grid_ComputePressure.o \
	Grid_ComputePressureDualEnergyFormalism.o \
	Grid_ComputeTemperatureField.o \
	Grid_ComputeTimeStep.o \
	Grid_IdentifyPhysicalQuantities.o \
	Grid_IdentifySpeciesFields.o \
	Grid_SetMinimumSupport.o \
	Grid_SolveHydroEquations.o \
	calcdiss.o \
	calc_dt.o \
	euler.o \
	f_message.o \
	inteuler.o \
	initialize_hydro.o \
	initialize_implosion.o \
	intvar.o \
	pgas2d_dual.o \
	pgas2d.o \
	ppm_de.o \
	sum_grid.o \
	twoshock.o \
	xeuler_sweep.o \
	yeuler_sweep.o \
	zeuler_sweep.o


FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX

FFLAGS   = $(FFLAGS_LOCAL)
CXXFLAGS = $(CXXFLAGS_LOCAL) $(INC_LOCAL)
LIB      = -lgfortran -lm


test_hydro: $(OBJS) test_hydro.o
	$(CXX) -o test_hydro test_hydro.o $(OBJS) $(LIB)

obj: $(OBJS)

test: test_hydro
	test_hydro

.PHONY: clean
clean:
	-@rm *.o test_hydro *~

include $(TOP)/src/MAKE/Make.include

include DEPEND
