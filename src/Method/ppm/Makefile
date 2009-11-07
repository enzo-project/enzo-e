TOP = ../../..

PACKAGE = method

TESTS = test_ppm

# This does no good--why not?  OMP_NUM_THREADS = 4 on stingmeyer
# PARALLEL = -fopenmp -ftree-parallelize-loops=4 

INC_LOCAL  = -I..
LIB_LOCAL  = $(LIB_HDF5) $(LIB_FORTRAN) $(PARALLEL)
CPPFLAGS_LOCAL   = -Dp8 -Dr8 -DSMALL_INTS -DLINUX $(PARALLEL)

OBJS_PPM = test_ppm.o \
           initialize_image.o \
           initialize_implosion.o \
           initialize_implosion3.o

OBJS = \
        calcdiss.o \
	calc_dt.o \
        cello_hydro.o \
        c_message.o \
	../CosmologyComputeExpansionFactor.o \
	../CosmologyComputeExpansionTimestep.o \
	CosmologyGetUnits.o \
	../data_dump.o \
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
	../Grid_SetExternalBoundaryValues.o \
	Grid_SetMinimumSupport.o \
	Grid_SolveHydroEquations.o \
	../initialize_hydro.o \
	inteuler.o \
	intvar.o \
	pgas2d_dual.o \
	pgas2d.o \
	ppm_de.o \
	../SetExternalBoundary.o \
	sum_field.o \
	twoshock.o \
	xeuler_sweep.o \
	yeuler_sweep.o \
	zeuler_sweep.o

new_all: includes all

include $(TOP)/src/MAKE/Make.include

test_ppm: includes $(OBJS) $(OBJS_PPM)
	$(CXX) -o test_ppm $(OBJS_PPM) $(OBJS) $(LIB)

initialize_image.o: image.h

image.h:
	ln -s circle2.h image.h

#-----------------------------------------------------------------------

.PHONY: includes
includes:
	-@ln -sf ../cello_hydro.h .
	-@ln -sf ../defines.h .

.PHONY: install-include
install-include:
	ln -sf `pwd`/cello_hydro.h             $(INCLUDE)

#-----------------------------------------------------------------------

.PHONY: clean-local
clean-local:
	-@rm -f *.png *.hdf5


include DEPEND
