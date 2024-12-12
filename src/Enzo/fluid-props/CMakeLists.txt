# See LICENSE_CELLO file for license and copyright information

# create Enzo::fluid-props
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_fluid-props
  fluid-props.hpp

  EnzoComputePressure.cpp EnzoComputePressure.hpp
  EnzoComputeTemperature.cpp EnzoComputeTemperature.hpp
  EnzoDualEnergyConfig.hpp
  EnzoEOSIdeal.hpp
  EnzoEOSIsothermal.hpp
  EnzoEOSVariant.cpp EnzoEOSVariant.hpp
  EnzoFluidFloorConfig.hpp
  EnzoPhysicsFluidProps.cpp EnzoPhysicsFluidProps.hpp
)

add_library(Enzo::fluid-props ALIAS Enzo_fluid-props)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_fluid-props PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_fluid-props PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_fluid-props PRIVATE ${Cello_TARGET_LINK_OPTIONS})
