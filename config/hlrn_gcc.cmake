# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)
if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for HLRN-IV supercomputer.\n"
    "This configuration has been tested using the following modules (last verified 2022-06-16):\n"
    "$ module load gcc/9.3.0 openmpi/gcc.9/3.1.5 hdf5/gcc.9/1.12.0 libpng/1.6.37 boost/1.72.0\n"
    "and with Charm++ configured as follows:\n"
    "$ cmake -DNETWORK=mpi -DSMP=OFF -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 ..\n"
    "Note: set PNG_ROOT to LIBPNG_PATH using -D or export.\n")

  set(CMAKE_CXX_COMPILER mpicxx CACHE STRING "")
  set(CMAKE_C_COMPILER mpicc CACHE STRING "")
  set(CMAKE_Fortran_COMPILER mpif90  CACHE STRING "")

  # the minimal set of required flags to successfully compile with this Fortran
  # compiler are handled internally (if those flags don't work, please update
  # the relevant internal logic rather than specifying them here)

  # add optional flags to C and C++ compilers that provide useful warnings
  #set(CMAKE_C_FLAGS "-Wall" CACHE STRING "Default C flags")
  #set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "Default C++ flags")

  # these flag(s) are currently only used when using openmp-simd optimizations
  # (to specify available/prefered instruction sets).
  # This particular value tells the compiler to optimize the code for the
  # instruction set of the machine used to compile the code.
  set(CONFIG_ARCH_FLAGS "-march=skylake-avx512")

  # Mark done
  set(__processedUserDefaults ON)

else()


endif()
