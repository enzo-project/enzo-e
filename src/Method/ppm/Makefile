TOP = ../../..

PACKAGE = method

TESTS = test_ppm test_implosion test_implosion3 test_image test_color

# This does no good--why not?  OMP_NUM_THREADS = 4 on stingmeyer
# PARALLEL = -fopenmp -ftree-parallelize-loops=4 

INC_LOCAL  = $(INC_HDF5)
LIB_LOCAL  = $(LIB_HDF5) $(LIB_FORTRAN) $(PARALLEL)
CPPFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX $(PARALLEL)

OBJS_PPM        = test_ppm.o initialize_ppm.o
OBJS_IMAGE      = test_image.o initialize_image.o
OBJS_COLOR      = test_color.o initialize_color.o
OBJS_IMPLOSION  = test_implosion.o  initialize_implosion.o
OBJS_IMPLOSION3 = test_implosion3.o initialize_implosion3.o

OBJS = \
        calcdiss.o \
	calc_dt.o \
        cello_hydro.o \
        c_message.o \
	CosmologyComputeExpansionFactor.o \
	CosmologyComputeExpansionTimestep.o \
	CosmologyGetUnits.o \
	data_dump.o \
	euler.o \
	FindField.o \
	f_message.o \
	Grid_ComputeGammaField.o \
	Grid_ComputePressureDualEnergyFormalism.o \
	Grid_ComputePressure.o \
	Grid_ComputeTemperatureField.o \
	Grid_ComputeTimeStep.o \
	Grid_IdentifyPhysicalQuantities.o \
	Grid_IdentifySpeciesFields.o \
	Grid_SetExternalBoundaryValues.o \
	Grid_SetMinimumSupport.o \
	Grid_SolveHydroEquations.o \
        initialize_hydro.o \
	inteuler.o \
	intvar.o \
	pgas2d_dual.o \
	pgas2d.o \
	ppm_de.o \
	print_field.o \
	SetExternalBoundary.o \
	sum_field.o \
	twoshock.o \
	xeuler_sweep.o \
	yeuler_sweep.o \
	zeuler_sweep.o

include $(TOP)/src/MAKE/Make.include

test_ppm: $(OBJS) $(OBJS_PPM)
	$(CXX) -o test_ppm $(OBJS_PPM) $(OBJS) $(LIB)
test_image: $(OBJS) $(OBJS_IMAGE)
	$(CXX) -o test_image $(OBJS_IMAGE) $(OBJS) $(LDFLAGS)
test_color: $(OBJS) $(OBJS_COLOR)
	$(CXX) -o test_color $(OBJS_COLOR) $(OBJS) $(LDFLAGS)
test_implosion: $(OBJS) $(OBJS_IMPLOSION)
	$(CXX) -o test_implosion $(OBJS_IMPLOSION) $(OBJS)  $(LDFLAGS)
test_implosion3: $(OBJS) $(OBJS_IMPLOSION3)
	$(CXX) -o test_implosion3 $(OBJS_IMPLOSION3) $(OBJS)  $(LDFLAGS)

.PHONY: clean-local
clean-local:
	-@rm -f *.png *.hdf5


include DEPEND
