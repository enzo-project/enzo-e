# See LICENSE_CELLO file for license and copyright information

# define the Enzo:hydro-mhd_toolkit subcomponent
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_hydro-mhd_toolkit
  EnzoSourceInternalEnergy.cpp EnzoSourceInternalEnergy.hpp
  EnzoSourceGravity.cpp EnzoSourceGravity.hpp
  EnzoLazyPassiveScalarFieldList.cpp EnzoLazyPassiveScalarFieldList.hpp
  EnzoBfieldMethod.cpp EnzoBfieldMethod.hpp
  EnzoBfieldMethodCT.cpp EnzoBfieldMethodCT.hpp
  EnzoIntegrationQuanUpdate.cpp
  EnzoReconstructor.cpp EnzoReconstructor.hpp
  EnzoReconstructorNN.cpp EnzoReconstructorNN.hpp
  EnzoReconstructorPLM.hpp
)

add_library(Enzo::hydro-mhd_toolkit ALIAS Enzo_hydro-mhd_toolkit)

# the following links a LOT of unnecessary stuff!
target_link_libraries(Enzo_hydro-mhd_toolkit PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_hydro-mhd_toolkit PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_hydro-mhd_toolkit PRIVATE ${Cello_TARGET_LINK_OPTIONS})
