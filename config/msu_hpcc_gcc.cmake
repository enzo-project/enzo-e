message(STATUS "Loading machine configuration for MSU's HPCC machine (incl. potential CUDA) builds.\n"
  "This configuration has been tested using the following modules (last verified 2021-07-14):\n"
  "$ module load -* gcccuda/2020a OpenMPI git Python HDF5 Boost libpng zlib/1.2.9\n")


set(CMAKE_CXX_COMPILER g++ CACHE STRING "")
set(CMAKE_C_COMPILER gcc CACHE STRING "")
set(CMAKE_Fortran_COMPILER gfortran CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-ffixed-line-length-132 -fdefault-real-8 -fdefault-double-8" CACHE STRING "Default Fortran flags")
