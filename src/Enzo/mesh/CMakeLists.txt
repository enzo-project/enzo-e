# See LICENSE_CELLO file for license and copyright information

# add source files related to amr to the enzo target. We explicitly list files
# since that make CMake more robust at rebuilds (after changing branches)
add_library(Enzo_mesh
  mesh.hpp
  EnzoProlong.cpp EnzoProlong.hpp
  EnzoProlongMC1.cpp EnzoProlongMC1.hpp
  EnzoProlongPoisson.cpp EnzoProlongPoisson.hpp
  EnzoRefineMass.cpp EnzoRefineMass.hpp
  EnzoRefineParticleMass.cpp EnzoRefineParticleMass.hpp
  EnzoRefineShock.cpp EnzoRefineShock.hpp
  EnzoRestrict.cpp EnzoRestrict.hpp
  interpolate.F
)
add_library(Enzo::mesh ALIAS Enzo_mesh)

# the following can be improved
target_link_libraries(Enzo_mesh PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_mesh PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_mesh PRIVATE ${Cello_TARGET_LINK_OPTIONS})

