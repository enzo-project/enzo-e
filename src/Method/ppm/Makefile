TOP = ../../..

TESTS = test_implosion test_implosion3 test_image test_implosion_opt

INC_LOCAL  = -I$(TOP)/include
LIB_LOCAL  = -L$(TOP)/lib -lcello

OBJS_IMAGE      = test_image.o initialize_image.o
OBJS_IMPLOSION  = test_implosion.o  initialize_implosion.o
OBJS_IMPLOSION3 = test_implosion3.o initialize_implosion3.o

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

# all: # WORKS ONLY IF empty all: TARGET IS HERE--WHY?

include $(TOP)/src/MAKE/Make.include

FFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -g $(DEF_OPTIMIZE)
CXXFLAGS_LOCAL = -Dp8 -Dr8 -DSMALL_INTS -DLINUX -g $(DEF_OPTIMIZE)

FFLAGS   = $(FFLAGS_LOCAL)
CXXFLAGS = $(CXXFLAGS_LOCAL) $(INC_LOCAL)
LIB      = $(LIB_LOCAL) -lgfortran -lm $(LIB_HDF5) $(PG)

test_image: $(OBJS) $(OBJS_FC) $(OBJS_IMAGE)
	$(CXX) -o test_image $(OBJS_IMAGE) $(OBJS) $(OBJS_FC) $(LIB)
test_implosion: $(OBJS) $(OBJS_FC) $(OBJS_IMPLOSION)
	$(CXX) -o test_implosion $(OBJS_IMPLOSION) $(OBJS) $(OBJS_FC)  $(LIB)
test_implosion3: $(OBJS) $(OBJS_FC) $(OBJS_IMPLOSION3)
	$(CXX) -o test_implosion3 $(OBJS_IMPLOSION3) $(OBJS) $(OBJS_FC)  $(LIB)
test_implosion_opt: $(OBJS) $(OBJS_FC_OPT) $(OBJS_IMPLOSION)
	$(CXX) -o test_implosion_opt $(OBJS_IMPLOSION) $(OBJS) $(OBJS_FC_OPT) $(LIB)

.PHONY: clean-local
clean-local:
	-@rm -f *.png *.hdf5


include DEPEND
