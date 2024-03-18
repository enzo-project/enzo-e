# See LICENSE_CELLO file for license and copyright information

# add  to the enzo target

# create Enzo::assorted, which is built from uncategorized source files.
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_assorted
  assorted.hpp
  EnzoMethodHeat.cpp EnzoMethodHeat.hpp
  EnzoMethodM1Closure.cpp EnzoMethodM1Closure.hpp
  EnzoMethodTurbulence.cpp EnzoMethodTurbulence.hpp
)
add_library(Enzo::assorted ALIAS Enzo_assorted)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_assorted PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_assorted PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_assorted PRIVATE ${Cello_TARGET_LINK_OPTIONS})
