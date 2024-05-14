# See LICENSE_CELLO file for license and copyright information

# create Enzo::cosmology
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_cosmology
  cosmology.hpp
  EnzoMethodComovingExpansion.hpp EnzoMethodComovingExpansion.cpp
  EnzoMethodCosmology.cpp EnzoMethodCosmology.hpp
  EnzoPhysicsCosmology.cpp EnzoPhysicsCosmology.hpp
  expand_terms.F
)

add_library(Enzo::cosmology ALIAS Enzo_cosmology)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_cosmology PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_cosmology PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_cosmology PRIVATE ${Cello_TARGET_LINK_OPTIONS})
