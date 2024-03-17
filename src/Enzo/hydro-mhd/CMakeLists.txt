# See LICENSE_CELLO file for license and copyright information

# STEP 1: define subcomponents
# (there's a lot in these subcomponents that the rest of Enzo-E doesn't need to
# know anything about)
add_subdirectory(ppm_fortran)  # define Enzo::ppm_fortran
add_subdirectory(ppml_fortran) # define Enzo::ppml_fortran
add_subdirectory(riemann)
add_subdirectory(toolkit) # define Enzo::hydro-mhd_toolkit

# STEP 2: define the Enzo::hydro-mhd target
#
# in the future, we may want to slightly refactor the files in this
# subdirectory so that they can be compiled into their own subtarget (there's
# a lot here that the rest of Enzo-E doesn't need to know anything about and
# this can improve compile times)

# Add source files related to hydro-mhd to the enzo target
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
  EnzoMHDIntegratorStageCommands.cpp EnzoMHDIntegratorStageCommands.hpp
  EnzoMethodPpm.cpp EnzoMethodPpm.hpp
)
add_library(Enzo::hydro-mhd ALIAS Enzo_hydro-mhd)

# the following can be improved
target_link_libraries(Enzo_hydro-mhd
  PUBLIC enzo ${CELLO_LIBS}
  PRIVATE Enzo::hydro-mhd_toolkit Enzo::ppm_fortran Enzo::ppml_fortran
          Enzo::riemann
)
target_include_directories(Enzo_hydro-mhd PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_hydro-mhd PRIVATE ${Cello_TARGET_LINK_OPTIONS})

