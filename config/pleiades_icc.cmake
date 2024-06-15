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

  # in principle we should set flags to specify hardware architecture (so that
  # they can be used with openmp-simd optimizations
  # set(CONFIG_ARCH_FLAGS ...)

  # if you choose to add other flags, you should generally prefer to use:
  #     ENZOE_C_FLIST_INIT, ENZOE_CXX_FLIST_INIT, ENZOE_Fortran_FLIST_INIT
  # rather than CMAKE_C_FLAGS, CMAKE_CXX_FLAGS, and CMAKE_Fortran_FLAGS
  # -> These alternatives will affect Cello/Enzo-E, but won't influence any
  #    dependencies compiled in the same-build
  # -> plus, the alternatives let users easily overwrite them

  # Set package paths (e.g., Grackle) - Only do this in personal machine files

  # Setting test environment
  set(PARALLEL_LAUNCHER "mpiexec" CACHE STRING "Use mpiexec for launching parallel tests")
  set(PARALLEL_LAUNCHER_NPROC_ARG "-np" CACHE STRING "Default option to set num procs")

  # Mark done
  set(__processedUserDefaults ON)

else()

endif()
