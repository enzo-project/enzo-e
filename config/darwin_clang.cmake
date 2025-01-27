# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)

# This is the configuration for a generic macOS installation using clang and
# gfortran
# - this was tested on a clean installation of version 12.5 of macOS on an
#   Intel-based MacBook Pro on July 24, 2022
# - historically there were challenges associated with installing Charm++ and
#   Cello/Enzo-E on Macs, but these issues seem to have been implicitly
#   resolved in each project's transition to using CMake for their
#   build-systems
#
# In this succesful installation:
# - Homebrew was used to install:
#       gfortran, hdf5, libpng, boost, and CMake
#   If you install this software through other channels, you may need to
#   provide CMake with hints about the locations of these packages (see the
#   main documentation for details).
# - the netlrts backend was used for Charm++
#
# If you're having problems with your installation:
# - did you make sure gfortran was installed before building Charm++?


if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for MacOS with clang and gfortran.\n"
    "If you encounter installation problems, please read the comments at the "
    "top of config/darwin_clang.cmake"
  )

  # Setting compilers
  set(CMAKE_CXX_COMPILER clang++ CACHE STRING "")
  set(CMAKE_C_COMPILER clang CACHE STRING "")
  set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")

  # the minimal set of required flags to successfully compile with this Fortran
  # compiler are handled internally (if those flags don't work, please update
  # the relevant internal logic rather than specifying them here)

  # these flag(s) are currently only used when using openmp-simd optimizations
  # (to specify available/prefered instruction sets).
  # This particular value tells the compiler to optimize the code for the
  # instruction set of the machine used to compile the code.
  set(CONFIG_ARCH_FLAGS "-march=native")

  # if you choose to add other flags, you should generally prefer to use:
  #     ENZOE_C_FLIST_INIT, ENZOE_CXX_FLIST_INIT, ENZOE_Fortran_FLIST_INIT
  # rather than CMAKE_C_FLAGS, CMAKE_CXX_FLAGS, and CMAKE_Fortran_FLAGS
  # -> These alternatives will affect Cello/Enzo-E, but won't influence any
  #    dependencies compiled in the same-build
  # -> plus, the alternatives let users easily overwrite them

  # in this case, we add a flag to unroll loops
  # (aside: this flag would be enabled anyways when OPTIMIZE_FP=TRUE)
  set(ENZOE_C_FLIST_INIT "-Wall;-funroll-loops")
  set(ENZOE_CXX_FLIST_INIT "${ENZOE_C_FLIST_INIT}")

  # Set package paths (e.g., Grackle) - Only do this in personal machine files

  # Mark done
  set(__processedUserDefaults ON)

else()

endif()
