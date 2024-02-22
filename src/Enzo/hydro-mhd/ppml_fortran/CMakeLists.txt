# See LICENSE_CELLO file for license and copyright information

# define the Enzo:ppml_fortran subcomponent
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_ppml_fortran
  ppml_fortran.hpp # included for completeness (but probably not necessary)

  PPML_Conservative.F
  PPML_HLLD.F
  PPML_MAIN.F
  PPML_monot.F
  PPML_potokz.F
  PPML_potokx.F
  PPML_potoky.F
  PPML_Primitive.F
  # ppml_fortran/PPML_QDD6.F # This is NOT used (can it be removed for good?)
  PPML_sub.F
  PPML_TimeStep.F
)
add_library(Enzo::ppml_fortran ALIAS Enzo_ppml_fortran)

# the following links a LOT of unnecessary stuff!
target_link_libraries(Enzo_ppml_fortran PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_ppml_fortran PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_ppml_fortran PRIVATE ${Cello_TARGET_LINK_OPTIONS})
