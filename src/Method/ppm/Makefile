TOP = ../../..

PACKAGE = method

TESTS = test_implosion test_implosion3 test_image test_color test_implosion_opt

# This does no good--why not?  OMP_NUM_THREADS = 4 on stingmeyer
# PARALLEL = -fopenmp -ftree-parallelize-loops=4 

INC_LOCAL  = $(INC_HDF5)
LIB_LOCAL  = $(LIB_HDF5) $(LIB_FORTRAN) $(PARALLEL)
CPPFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX $(DEF_OPTIMIZE) $(PARALLEL)

OBJS_IMAGE      = test_image.o initialize_image.o
OBJS_COLOR      = test_color.o initialize_color.o
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

include $(TOP)/src/MAKE/Make.include

test_image: $(OBJS) $(OBJS_FC) $(OBJS_IMAGE)
	$(CXX) -o test_image $(OBJS_IMAGE) $(OBJS) $(OBJS_FC) $(LDFLAGS)
test_color: $(OBJS) $(OBJS_FC) $(OBJS_COLOR)
	$(CXX) -o test_color $(OBJS_COLOR) $(OBJS) $(OBJS_FC) $(LDFLAGS)
test_implosion: $(OBJS) $(OBJS_FC) $(OBJS_IMPLOSION)
	$(CXX) -o test_implosion $(OBJS_IMPLOSION) $(OBJS) $(OBJS_FC)  $(LDFLAGS)
test_implosion3: $(OBJS) $(OBJS_FC) $(OBJS_IMPLOSION3)
	$(CXX) -o test_implosion3 $(OBJS_IMPLOSION3) $(OBJS) $(OBJS_FC)  $(LDFLAGS)
test_implosion_opt: $(OBJS) $(OBJS_FC_OPT) $(OBJS_IMPLOSION)
	$(CXX) -o test_implosion_opt $(OBJS_IMPLOSION) $(OBJS) $(OBJS_FC_OPT) $(LDFLAGS)

.PHONY: clean-local
clean-local:
	-@rm -f *.png *.hdf5


include DEPEND
