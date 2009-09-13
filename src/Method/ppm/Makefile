#
# test         Run unit tests
# clean        Clean directory
# test_array   Compile the unit tests
# dep          Generate dependencies
# 

TOP = ../../..

TEST_EXES = test_implosion test_enzo test_sdsc

INC_LOCAL  = -I$(TOP)/include
LIB_LOCAL  = -L$(TOP)/lib -lcello

OBJS_ENZO = test_enzo.o initialize_enzo.o
OBJS_SDSC = test_sdsc.o initialize_sdsc.o
OBJS_IMPLOSION = test_implosion.o initialize_implosion.o

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
	Grid_SetExternalBoundaryValues.o \
	SetExternalBoundary.o \
	calcdiss.o \
	calc_dt.o \
	euler.o \
	f_message.o \
	inteuler.o \
	initialize_hydro.o \
	intvar.o \
	pgas2d_dual.o \
	pgas2d.o \
	ppm_de.o \
	twoshock.o \
	xeuler_sweep.o \
	yeuler_sweep.o \
	zeuler_sweep.o \
	print_field.o \
	sum_field.o \
	data_dump.o


FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3

FFLAGS   = $(FFLAGS_LOCAL)
CXXFLAGS = $(CXXFLAGS_LOCAL) $(INC_LOCAL)
LIB      = $(LIB_LOCAL) -lgfortran -lm -lhdf5

all: test_enzo test_implosion test_sdsc
test_enzo: $(OBJS) $(OBJS_ENZO)
	$(CXX) -o test_enzo $(OBJS_ENZO) $(OBJS) $(LIB)
test_sdsc: $(OBJS) $(OBJS_SDSC) initialize_sdsc.o
	$(CXX) -o test_sdsc $(OBJS_SDSC) $(OBJS) $(LIB)
test_implosion: $(OBJS) $(OBJS_IMPLOSION)
	$(CXX) -o test_implosion $(OBJS_IMPLOSION) $(OBJS) $(LIB)

.PHONY: install
install: $(TEST_EXES)
	cp $(TEST_EXES) $(TOP)/bin

obj: $(OBJS)

test: $(TEST_EXES)
#	./test_implosion # slow
#	./test_enzo      # slow

.PHONY: clean
clean:
	-@rm *.o $(TEST_EXES) *~

include $(TOP)/src/MAKE/Make.include

include DEPEND
