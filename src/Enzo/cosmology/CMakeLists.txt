# See LICENSE_CELLO file for license and copyright information

# add source files related to cosmology to the enzo target

# Get the list of source files in this directory
# - we do this using GLOB patterns. This approach is not recommended by the
#   authors of CMake (their recommendation is to explicitly list all files that
#   must be installed).
# - Some of the disadvantages of this approach are mitigated by inclusion of
#   the CONFIGURE_DEPENDS flag.
# - See the CMake Primer section of the developer documentation for more details
file(GLOB LOCAL_SRC_FILES CONFIGURE_DEPENDS
  *.cpp *.hpp *F
)

target_sources(enzo PRIVATE ${LOCAL_SRC_FILES})
