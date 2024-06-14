# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)
if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for NASA Pleiades supercomputer.\n"
    "This configuration has been tested using the following modules:\n"
    "$ module use -a /nasa/modulefiles/testing\n"
    "$ module purge\n"
    "$ module load pkgsrc/2021Q1 gcc/10.2 mpi-hpe/mpt.2.23 boost/1.76 comp-intel/2020.4.304 hdf5/1.12.0_serial pkgsrc/2021Q1\n")

  set(CMAKE_CXX_COMPILER icpc CACHE STRING "")
  set(CMAKE_C_COMPILER icc CACHE STRING "")
  set(CMAKE_Fortran_COMPILER ifort CACHE STRING "")

  # the minimal set of required flags to successfully compile with this Fortran
  # compiler are handled internally (if those flags don't work, please update
  # the relevant internal logic rather than specifying them here)

  # Setting package paths (e.g., Grackle)

  # Setting test environment
  set(PARALLEL_LAUNCHER "mpiexec" CACHE STRING "Use mpiexec for launching parallel tests")
  set(PARALLEL_LAUNCHER_NPROC_ARG "-np" CACHE STRING "Default option to set num procs")

  # Mark done
  set(__processedUserDefaults ON)

else()

endif()
