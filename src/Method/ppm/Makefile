#
# test         Run unit tests
# clean        Clean directory
# test_array   Compile the unit tests
# dep          Generate dependencies
# 

# DEF_OPTIMIZE = -D OPTIMIZE_CELL_WIDTH
TOP = ../../..

PG = 

TEST_EXES = test_implosion test_image

INC_LOCAL  = -I$(TOP)/include
LIB_LOCAL  = -L$(TOP)/lib -lcello

OBJS_IMAGE     = test_image.o initialize_image.o
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
	print_field.o \
	sum_field.o \
	data_dump.o \
        initialize_hydro.o

OBJS_FC = calcdiss.o \
	calc_dt.o \
	euler.o \
	f_message.o \
	inteuler.o \
	intvar.o \
	pgas2d_dual.o \
	pgas2d.o \
	ppm_de.o \
	twoshock.o \
	xeuler_sweep.o \
	yeuler_sweep.o \
	zeuler_sweep.o

OBJS_FC_OPT = ppm.o

FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3 $(PG) $(DEF_OPTIMIZE)
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -O3 $(PG) $(DEF_OPTIMIZE)

FFLAGS   = $(FFLAGS_LOCAL)
CXXFLAGS = $(CXXFLAGS_LOCAL) $(INC_LOCAL)
LIB      = $(LIB_LOCAL) -lgfortran -lm -lhdf5 $(PG)

all: test_image test_implosion
test_image: $(OBJS) $(OBJS_FC) $(OBJS_IMAGE)
	$(CXX) -o test_image $(OBJS_IMAGE) $(OBJS) $(OBJS_FC) $(LIB)
test_implosion: $(OBJS) $(OBJS_FC) $(OBJS_IMPLOSION)
	$(CXX) -o test_implosion $(OBJS_IMPLOSION) $(OBJS) $(OBJS_FC)  $(LIB)
test_implosion_opt: $(OBJS) $(OBJS_FC_OPT) $(OBJS_IMPLOSION)
	$(CXX) -o test_implosion_opt $(OBJS_IMPLOSION) $(OBJS) $(OBJS_FC_OPT)  $(LIB)

.PHONY: install
install: $(TEST_EXES)
	cp $(TEST_EXES) $(TOP)/bin

obj: $(OBJS)

test:

.PHONY: clean
clean:
	-@rm *.o $(TEST_EXES) *~

include $(TOP)/src/MAKE/Make.include

include DEPEND
