# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)
if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for generic Linux machine with GCC stack.\n")

  # Setting compilers
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
  set(CMAKE_C_COMPILER gcc CACHE STRING "")
  set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")

  # the minimal set of required flags to successfully compile with this Fortran
  # compiler are handled internally (if those flags don't work, please update
  # the relevant internal logic rather than specifying them here)

  # Note (12/2021): passing -march=native to gfortran seems to slow down the
  # PPM solver

  set(__ARCH_C_OPT_FLAGS "-O0 -g -funroll-loops")
  set(CMAKE_C_FLAGS_RELEASE "${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${__ARCH_C_OPT_FLAGS}")

  # these flag(s) are currently only used when using openmp-simd optimizations
  # (to specify available/prefered instruction sets).
  # This particular value tells the compiler to optimize the code for the
  # instruction set of the machine used to compile the code.
  set(CONFIG_ARCH_FLAGS "-march=native -g")

  # if you choose to add other flags, you should generally prefer to use:
  #     ENZOE_C_FLIST_INIT, ENZOE_CXX_FLIST_INIT, ENZOE_Fortran_FLIST_INIT
  # rather than CMAKE_C_FLAGS, CMAKE_CXX_FLAGS, and CMAKE_Fortran_FLAGS
  # -> These alternatives will affect Cello/Enzo-E, but won't influence any
  #    dependencies compiled in the same-build
  # -> plus, the alternatives let users easily overwrite them

  # add optional flags to C and C++ compilers that provide useful warnings
  # (also add flag to unroll loops - this flag would also be enabled anyways
  # when OPTIMIZE_FP=TRUE)
  set(ENZOE_C_FLIST_INIT "-Wall;-funroll-loops")
  set(ENZOE_CXX_FLIST_INIT "${ENZOE_C_FLIST}")

  # Set package paths (e.g., Grackle) - Only do this in personal machine files

  # Mark done
  set(__processedUserDefaults ON)

else()


endif()
