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

  # the minimal set of required flags to successfully compile with this Fortran
  # compiler are handled internally (if those flags don't work, please update
  # the relevant internal logic rather than specifying them here)

  # these flag(s) are currently only used when using openmp-simd optimizations
  # (to specify available/prefered instruction sets).
  # This particular value tells the compiler to optimize the code for the
  # instruction set of the machine used to compile the code.
  set(CONFIG_ARCH_FLAGS "-xHost")

  # if you choose to add other flags, you should generally prefer to use:
  #     ENZOE_C_FLIST_INIT, ENZOE_CXX_FLIST_INIT, ENZOE_Fortran_FLIST_INIT
  # rather than CMAKE_C_FLAGS, CMAKE_CXX_FLAGS, and CMAKE_Fortran_FLAGS
  # -> These alternatives will affect Cello/Enzo-E, but won't influence any
  #    dependencies compiled in the same-build
  # -> plus, the alternatives let users easily overwrite them

  # add optional flags to C and C++ compilers that provide useful warnings
  set(ENZOE_C_FLIST_INIT "-Wall")
  set(ENZOE_CXX_FLIST_INIT "${ENZOE_C_FLIST_INIT}")

  # Set package paths (e.g., Grackle) - Only do this in personal machine files

  # Mark done
  set(__processedUserDefaults ON)

else()

endif()
