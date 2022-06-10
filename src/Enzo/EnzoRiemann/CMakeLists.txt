# The riemann library implements the C++ Riemann Solvers
#
# The main purpose of this library is to isolate the template-heavy code used
# implement the Riemann Solvers from the rest of the Enzo layer. Only changes
# to the public header file, EnzoRiemann.hpp, should require the rest of the
# Enzo layer to be rebuilt

file(GLOB RIEMANN_SRC_FILES ./EnzoRiemann*.cpp ./EnzoRiemann*.hpp)
add_library(riemann STATIC ${RIEMANN_SRC_FILES})
# the following line is *ONLY* required because this library includes
#   src/Enzo/enzo.hpp which effectively includes all of the other header files
#   in the Enzo/Cello layers. Since some of these may include headers from
#   external libraries (e.g. src/Cello/_disk.hpp includes <hdf5.h>), we need to
#   tell the CMake to link against these external libraries so that it can
#   inform the compiler to search for headers in appropriate directories for
#   each dependency 
target_link_libraries(riemann PRIVATE ${External_LIBS})

add_dependencies(riemann
  enzoCharmModule
  main_enzoCharmModule
  simulationCharmModule
)
# this could probably be paired down
target_include_directories (riemann
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} # this is public to inform the compiler where to
                              # search for the public header file,
                              # "EnzoRiemann.hpp"
                              # (also specifies the location of the library's
                              # private header files)
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/.. # for including "enzo.hpp"
  ${CMAKE_CURRENT_BINARY_DIR}/.. # for including "enzo.decl.hpp" (implicitly
                                 # included by "enzo.hpp")
  ${CMAKE_CURRENT_SOURCE_DIR}/../../Cello # for including "cello.hpp"
  ${CMAKE_CURRENT_BINARY_DIR}/../../Cello # for including "simulation.decl.h"
                                          # (implicitly included by "cello.hpp")
  ${CHARM_INCLUDE_DIRS} # for including <charm++.h> (implicitly included by
                        # "cello.hpp")
)
target_link_options(riemann PRIVATE ${Cello_TARGET_LINK_OPTIONS})
