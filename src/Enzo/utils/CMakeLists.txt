# See LICENSE_CELLO file for license and copyright information

# create Enzo::tests, which represents defines source files for initializers
# that explicitly exist for testing purposes
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_utils
  utils.hpp utils.cpp
  EnzoCenteredFieldRegistry.cpp EnzoCenteredFieldRegistry.hpp
  EnzoComputeCicInterp.cpp EnzoComputeCicInterp.hpp
  EnzoEFltArrayMap.cpp EnzoEFltArrayMap.hpp
  EnzoFieldAdaptor.cpp EnzoFieldAdaptor.hpp
)
add_library(Enzo::utils ALIAS Enzo_utils)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_utils PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_utils PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_utils PRIVATE ${Cello_TARGET_LINK_OPTIONS})
