# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)
if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for generic Linux machine with GCC stack.\n")

  # Setting compilers
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
  set(CMAKE_C_COMPILER gcc CACHE STRING "")
  set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")
  set(CMAKE_Fortran_FLAGS "-ffixed-line-length-132" CACHE STRING "Default Fortran flags")

  # Setting package paths (e.g., Grackle)

  # Mark done
  set(__processedUserDefaults ON)

else()

  if (USE_DOUBLE_PREC)
    string(APPEND CMAKE_Fortran_FLAGS " -fdefault-real-8 -fdefault-double-8")
  endif()

endif()
