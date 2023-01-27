# See LICENSE_CELLO file for license and copyright information

# STEP 1: adds source files related to Hydro/MHD integrators to the enzo target
#
# in the future, we may want to slightly refactor the files in this
# subdirectory so that they can be compiled into their own subtarget (there's
# a lot here that the rest of Enzo-E doesn't need to know anything about and
# this can improve compile times)

# Get the list of source files in this directory & the fortran subdirectories
# - we do this using GLOB patterns. This approach is not recommended by the
#   authors of CMake (their recommendation is to explicitly list all files that
#   must be installed).
# - Some of the disadvantages of this approach are mitigated by inclusion of
#   the CONFIGURE_DEPENDS flag.
# - See the CMake Primer section of the developer documentation for more details
file(GLOB LOCAL_SRC_FILES CONFIGURE_DEPENDS
  *.cpp *.hpp ppm_fortran/*.F ppml_fortran/*.F
)

# TODO: Check if that file can be removed for good. It also not being used in the SCons build.
# -> PPML_sub.F:(.text+0x5ebf): multiple definition of `qdd6_'; CMakeFiles/enzo-e.dir/PPML_QDD6.F.o:PPML_QDD6.F:(.text+0x0): first defined here
list(FILTER LOCAL_SRC_FILES EXCLUDE REGEX "PPML_QDD6")

target_sources(enzo PRIVATE ${LOCAL_SRC_FILES})

# STEP 2: define the riemann library
add_subdirectory(riemann)
