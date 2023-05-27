# See LICENSE_CELLO file for license and copyright information

# add source files related to heating and cooling to the enzo target
#
# in the near future, it's my intention to make most of the contents in this
# directory into their own library and make it the only part of Enzo-E that
# depends on Grackle (this will make incremental rebuilds that turn Grackle
# on/off a lot faster)

# Get the list of source files in this directory
# - we do this using GLOB patterns. This approach is not recommended by the
#   authors of CMake (their recommendation is to explicitly list all files that
#   must be installed).
# - Some of the disadvantages of this approach are mitigated by inclusion of
#   the CONFIGURE_DEPENDS flag.
# - See the CMake Primer section of the developer documentation for more details
file(GLOB LOCAL_SRC_FILES CONFIGURE_DEPENDS
  *.cpp *.hpp
)

target_sources(enzo PRIVATE ${LOCAL_SRC_FILES})
