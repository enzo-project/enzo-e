# See LICENSE_CELLO file for license and copyright information

# add source files related to the enzo target
#
# NOTE: this subdirectory was originally just named core, but that caused
#       charmc to produce an error message when it was called to translate a
#       charm++ interface file.
#       - The error message is benign. It occurs during cleanup in charmc when
#         it calls `rm -f core` (presumably to remove any core dump files), and
#         rm complains that it can't remove a directory.
#       - We choose to avoid that name in order to avoid worrying users

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

# remove the unit-test file from this search
list(FILTER LOCAL_SRC_FILES EXCLUDE REGEX "test_EnzoUnits")

target_sources(enzo PRIVATE ${LOCAL_SRC_FILES})

# Add a unit test
add_executable(test_enzo_units "test_EnzoUnits.cpp")
target_link_libraries(test_enzo_units PRIVATE enzo main_enzo)
target_link_options(test_enzo_units PRIVATE ${Cello_TARGET_LINK_OPTIONS})
