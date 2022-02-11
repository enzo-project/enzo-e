# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)
if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for linux icc.\n"
    "This configuration has been tested (2021-11-13) using the intel compilers\n"
    "downloaded as part of the Intel oneAPI HPC Toolkit (the particular toolkit\n"
    "was from January 2021).\n"
    "For this build, Charm++ was configured as follows:\n"
    "$ cmake -DNETWORK=netlrts -DSMP=OFF -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort ..\n")

  # first, set the compilers and mandatory flags
  set(CMAKE_CXX_COMPILER icpc CACHE STRING "")
  set(CMAKE_C_COMPILER icc CACHE STRING "")
  set(CMAKE_Fortran_COMPILER ifort CACHE STRING "")
  set(CMAKE_Fortran_FLAGS "-nofor-main" CACHE STRING "Default Fortran flags")

  # add optional flags to C and C++ compilers that provide useful warnings
  set(CMAKE_C_FLAGS "-Wall" CACHE STRING "Default C flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "Default C++ flags")

  # these flag(s) are currently only used when using openmp-simd optimizations
  # (to specify available/prefered instruction sets).
  # This particular value tells the compiler to optimize the code for the
  # instruction set of the machine used to compile the code.
  set(CONFIG_ARCH_FLAGS "-xHost")

  # Setting package paths (e.g., Grackle)

  # Mark done
  set(__processedUserDefaults ON)

else()

  if (USE_DOUBLE_PREC)
    string(APPEND CMAKE_Fortran_FLAGS " -real-size 64 -double-size 64")
  endif()

endif()
