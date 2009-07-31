#
# test         Run unit tests
# clean        Clean directory
# test_array   Compile the unit tests
# dep          Generate dependencies
# 

TOP = ../../..

include $(TOP)/src/MAKE/Make.include

INC_LOCAL  = -I$(TOP)/include

OBJS = \
        archsinh.o \
	calcdiss.o \
	CosmologyComputeExpansionFactor.o \
	CosmologyGetUnits.o \
	CosmologyParameters.o \
	euler.o \
        c_message.o \
        f_message.o \
        global_data.o \
	FindField.o \
        Gadget.o \
	Grid_ComputeGammaField.o \
	Grid_ComputePressure.o \
	Grid_ComputePressureDualEnergyFormalism.o \
	Grid_ComputeTemperatureField.o \
	Grid_DebugCheck.o \
	Grid_GadgetComputeTemperature.o \
	Grid_GadgetComputeTemperatureDEF.o \
	Grid_Gadgetconvertutotemp.o \
	Grid_Gadgetfindabundancesandrates.o \
	Grid_IdentifyPhysicalQuantities.o \
	Grid_IdentifySpeciesFields.o \
	Grid_SetMinimumSupport.o \
	Grid_SolveHydroEquations.o \
	inteuler.o \
	intvar.o \
	pgas2d_dual.o \
	pgas2d.o \
	ppm_de.o \
	twoshock.o \
	xeuler_sweep.o \
	yeuler_sweep.o \
	zeuler_sweep.o

FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX

FFLAGS   = $(FFLAGS_LOCAL)
CXXFLAGS = $(CXXFLAGS_LOCAL) $(INC_LOCAL)
LIB      = -lgfortran -lm


obj: $(OBJS)

test_hydro: $(OBJS) test_hydro.o
	$(CXX) -o test_hydro test_hydro.o $(OBJS) $(LIB)

test: test_hydro
	test_hydro

.PHONY: clean
clean:
	-rm *.o test_hydro *~ out.make.dep Makefile.bak

