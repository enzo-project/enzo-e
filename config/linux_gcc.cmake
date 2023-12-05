# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)
if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for generic Linux machine with GCC stack.\n")

  # Remove cmake -DNDEBUG -O3
  set(CMAKE_C_FLAGS_INIT "")
  set(CMAKE_CXX_FLAGS_INIT "")
  set(CMAKE_Fortran_FLAGS_INIT "")

  # Setting compilers
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
  set(CMAKE_C_COMPILER gcc CACHE STRING "")
  set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")
  # Note (12/2021): passing -march=native to gfortran seems to slow down the
  # PPM solver
  set(CMAKE_Fortran_FLAGS "-ffixed-line-length-132" CACHE STRING "Default Fortran flags")

  set(__ARCH_C_OPT_FLAGS "-O3 -g -funroll-loops")
  set(CMAKE_C_FLAGS_RELEASE "${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${__ARCH_C_OPT_FLAGS}")

  # these flag(s) are currently only used when using openmp-simd optimizations
  # (to specify available/prefered instruction sets).
  # This particular value tells the compiler to optimize the code for the
  # instruction set of the machine used to compile the code.
  set(CONFIG_ARCH_FLAGS "-march=native")

  # add optional flags to C and C++ compilers that provide useful warnings
  # (also add flag to unroll loops - this flag would also be enabled anyways
  # when OPTIMIZE_FP=TRUE)
  set(CMAKE_C_FLAGS "-Wall -funroll-loops" CACHE STRING "Default C flags")
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "Default C++ flags")

  # Setting package paths (e.g., Grackle) - (meant for personal machine files)

  # Mark done
  set(__processedUserDefaults ON)

else()

  if (USE_DOUBLE_PREC)
    string(APPEND CMAKE_Fortran_FLAGS " -fdefault-real-8 -fdefault-double-8")
  endif()

endif()
