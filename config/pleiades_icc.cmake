message(STATUS "Loading machine configuration for NASA Pleiades supercomputer.\n"
  "This configuration has been tested using the following modules:\n"
  "$ module use -a /nasa/modulefiles/testing\n"
  "$ module purge\n"
  "$ module load pkgsrc/2021Q1 gcc/10.2 mpi-hpe/mpt.2.23 boost/1.75 comp-intel/2020.4.304 hdf5/1.12.0_serial pkgsrc/2021Q1\n")


set(CMAKE_CXX_COMPILER icpc CACHE STRING "")
set(CMAKE_C_COMPILER icc CACHE STRING "")
set(CMAKE_Fortran_COMPILER ifort CACHE STRING "")
set(CMAKE_Fortran_FLAGS "-real-size 64 -double-size 64 -nofor-main" CACHE STRING "Default Fortran flags")
