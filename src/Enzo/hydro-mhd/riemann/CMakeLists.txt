# The riemann library implements the C++ Riemann Solvers
#
# The main purpose of this library is to isolate the template-heavy code used
# to implement the Riemann Solvers from the rest of the Enzo layer. Under the
# current setup:
# - code outside this directory will only need to be recompiled when the public
#   header file, EnzoRiemann.hpp, is mutated.
# - Changes to any other file will only trigger re-linking

# Get the list of source files
# - we do this using GLOB patterns. This approach is not recommended by the
#   authors of CMake (their recommendation is to explicitly list all files that
#   must be installed).
# - Some of the disadvantages of this approach are mitigated by inclusion of
#   the CONFIGURE_DEPENDS flag.
# - See the CMake Primer section of the developer documentation for more detail
file(GLOB RIEMANN_SRC_FILES CONFIGURE_DEPENDS
  ./EnzoRiemann*.cpp ./EnzoRiemann*.hpp
)

add_library(riemann STATIC ${RIEMANN_SRC_FILES})

# the following line is *ONLY* required because the source file includes
# "cello.hpp" and "enzo.hpp".
# - now that things are set up appropriately, the include-directories for
#   transitive dependencies (inherited via headers accessed through "cello.hpp"
#   or "enzo.hpp") are properly inherited.
# - in the future, it might be nice to only depend on subcomponents of the
#   Cello or Enzo layer (to reduce how often this library gets recompiled)
target_link_libraries(riemann PRIVATE enzo ${Cello_LIBS})

# the include_directory is public to inform dependents where to search for the
# public header file, "EnzoRiemann.hpp".
# (It's also used during compilation of the library since it also specifies the
# directory holding the private header files)
target_include_directories (riemann PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_options(riemann PRIVATE ${Cello_TARGET_LINK_OPTIONS})
