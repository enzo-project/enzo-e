# See LICENSE_CELLO file for license and copyright information

# add source files for code involving enzo-specific io to the enzo target

# create Enzo::io, which includes source files enzo-specific io
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_io
  io.hpp

  control_restart.cpp
  EnzoInitialHdf5.cpp
  EnzoInitialHdf5.hpp
  EnzoInitialMusic.cpp
  EnzoInitialMusic.hpp
  EnzoMethodCheck.cpp
  EnzoMethodCheck.hpp
  IoEnzoBlock.cpp
  IoEnzoBlock.hpp
  IoEnzoReader.hpp
  IoEnzoWriter.hpp
)

add_library(Enzo::io ALIAS Enzo_io)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_io PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_io PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_io PRIVATE ${Cello_TARGET_LINK_OPTIONS})
