# See LICENSE_CELLO file for license and copyright information

# create Enzo::chemistry
#
# NOTE: in the long-term, the goal is to isolate Grackle as a private
#       dependency of this target. By doing that, we can rebuild Enzo
#       with/without linking against Grackle without needing to recompile any
#       source files other than the ones used in this target.
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_chemistry
  chemistry.hpp
  EnzoComputeCoolingTime.cpp EnzoComputeCoolingTime.hpp
  EnzoMethodGrackle.cpp EnzoMethodGrackle.hpp
  GrackleChemistryData.cpp GrackleChemistryData.hpp
  GrackleFacade.cpp GrackleFacade.hpp
)

add_library(Enzo::chemistry ALIAS Enzo_chemistry)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_chemistry PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_chemistry PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_chemistry PRIVATE ${Cello_TARGET_LINK_OPTIONS})
