# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)
if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for generic Linux machine with GCC stack.\n")

  # Setting compilers
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
  set(CMAKE_C_COMPILER gcc CACHE STRING "")
  set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")
  # Note (12/2021): passing -march=native to gfortran seems to slow down the
  # PPM solver
  set(CMAKE_Fortran_FLAGS "-ffixed-line-length-132" CACHE STRING "Default Fortran flags")

  # Set some architecture-specific optimization flags
  # (in the future, we might want to excercise more control over the targetted
  # instruction set to a higher level)

  # if we use -fno-math-errno, that cause 1D MHD shock tube problems to have
  # slightly different L1 error norms for the VL+CT solver, depending on the
  # position in the transverse direction (see MHD_shock_tube_test)
  set(__ARCH_C_OPT_FLAGS "-O3 -fopenmp-simd -march=native -DNDEBUG -funroll-loops -fno-trapping-math -fno-signed-zeros")

  set(CMAKE_C_FLAGS_RELEASE "${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g ${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${__ARCH_C_OPT_FLAGS}")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g ${__ARCH_C_OPT_FLAGS}")

  # Setting package paths (e.g., Grackle)

  # Mark done
  set(__processedUserDefaults ON)

else()

  if (USE_DOUBLE_PREC)
    string(APPEND CMAKE_Fortran_FLAGS " -fdefault-real-8 -fdefault-double-8")
  endif()

endif()
