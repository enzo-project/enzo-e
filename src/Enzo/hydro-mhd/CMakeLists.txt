# See LICENSE_CELLO file for license and copyright information

# STEP 1: define subcomponents
# (there's a lot in these subcomponents that the rest of Enzo-E doesn't need to
# know anything about)
add_subdirectory(riemann)
# in the future, we may probably want to define each subdirectory as a
# subcomponent

# STEP 2: define the Enzo::hydro-mhd target
#
# in the future, we may want to slightly refactor the files in this
# subdirectory so that they can be compiled into their own subtarget (there's
# a lot here that the rest of Enzo-E doesn't need to know anything about and
# this can improve compile times)

# Add source files related to problem-initialization to the enzo target
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_hydro-mhd

  hydro-mhd.hpp
  EnzoMethodMHDVlct.cpp EnzoMethodMHDVlct.hpp
  EnzoMethodPpml.cpp EnzoMethodPpml.hpp
  SetMinimumSupport.cpp
  SolveMHDEquations.cpp
  SolveHydroEquations.cpp
  EnzoMethodPpm.cpp EnzoMethodPpm.hpp

  toolkit/EnzoSourceInternalEnergy.cpp toolkit/EnzoSourceInternalEnergy.hpp
  toolkit/EnzoSourceGravity.cpp toolkit/EnzoSourceGravity.hpp
  toolkit/EnzoLazyPassiveScalarFieldList.cpp toolkit/EnzoLazyPassiveScalarFieldList.hpp
  toolkit/EnzoBfieldMethod.cpp toolkit/EnzoBfieldMethod.hpp
  toolkit/EnzoBfieldMethodCT.cpp toolkit/EnzoBfieldMethodCT.hpp
  toolkit/EnzoIntegrationQuanUpdate.cpp
  toolkit/EnzoReconstructor.cpp toolkit/EnzoReconstructor.hpp
  toolkit/EnzoReconstructorNN.cpp toolkit/EnzoReconstructorNN.hpp
  toolkit/EnzoReconstructorPLM.hpp

  ppm_fortran/flux_hllc.F
  ppm_fortran/interp1d.F
  ppm_fortran/int_lin3d.F
  ppm_fortran/cicinterp.F
  ppm_fortran/twoshock.F
  ppm_fortran/ngpinterp.F
  ppm_fortran/flux_hll.F
  ppm_fortran/intprim.F
  ppm_fortran/calcdiss.F
  ppm_fortran/tscint2d.F
  ppm_fortran/ppm_de.F
  ppm_fortran/cic_deposit.F
  ppm_fortran/flux_twoshock.F
  ppm_fortran/interp3d.F
  ppm_fortran/pgas2d.F
  ppm_fortran/cic_interp.F
  ppm_fortran/euler.F
  ppm_fortran/intpos.F
  ppm_fortran/inteuler.F
  ppm_fortran/yeuler_sweep.F
  ppm_fortran/calc_eigen.F
  ppm_fortran/pgas2d_dual.F
  ppm_fortran/xeuler_sweep.F
  ppm_fortran/calc_dt.F
  ppm_fortran/intvar.F
  ppm_fortran/tscint3d.F
  ppm_fortran/zeuler_sweep.F
  ppm_fortran/interp2d.F
  ppm_fortran/tscint1d.F

  ppml_fortran/PPML_potoky.F
  ppml_fortran/PPML_potokx.F
  ppml_fortran/PPML_HLLD.F
  ppml_fortran/PPML_monot.F
  # ppml_fortran/PPML_QDD6.F # This is NOT used (can it be removed for good?)
  ppml_fortran/PPML_potokz.F
  ppml_fortran/PPML_MAIN.F
  ppml_fortran/PPML_Primitive.F
  ppml_fortran/PPML_sub.F
  ppml_fortran/PPML_TimeStep.F
  ppml_fortran/PPML_Conservative.F
)
add_library(Enzo::hydro-mhd ALIAS Enzo_hydro-mhd)

# the following can be improved
target_link_libraries(Enzo_hydro-mhd
  PUBLIC enzo ${CELLO_LIBS}
  PRIVATE riemann
)
target_include_directories(Enzo_hydro-mhd PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_hydro-mhd PRIVATE ${Cello_TARGET_LINK_OPTIONS})

