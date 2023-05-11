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
  set(CMAKE_Fortran_FLAGS "-ffixed-line-length-132" CACHE STRING "Default Fortran flags")

  # these flag(s) are currently only used when using openmp-simd optimizations
  # (to specify available/prefered instruction sets).
  # This particular value tells the compiler to optimize the code for the
  # instruction set of the machine used to compile the code.
  set(CONFIG_ARCH_FLAGS "-march=native")

  # add flag to unroll loops (this flag would also be enabled anyways when
  # OPTIMIZE_FP=TRUE)
  set(CMAKE_C_FLAGS "-funroll-loops" CACHE STRING "Default C flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "Default C++ flags")

  # Setting package paths (e.g., Grackle) - (meant for personal machine files)

  # Mark done
  set(__processedUserDefaults ON)

else()

  if (USE_DOUBLE_PREC)
    string(APPEND CMAKE_Fortran_FLAGS " -fdefault-real-8 -fdefault-double-8")
  endif()

endif()
