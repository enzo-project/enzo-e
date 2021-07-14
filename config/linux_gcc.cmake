message(STATUS "Loading machine configuration for generic Linux machine with GCC stack.\n")

set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
set(CMAKE_C_COMPILER gcc CACHE STRING "")
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-ffixed-line-length-132 -fdefault-real-8 -fdefault-double-8" CACHE STRING "Default Fortran flags")
