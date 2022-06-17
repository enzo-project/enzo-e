# Machine configs are included twice.
# First to set compilers, paths, default options.
# Second to set dependent options (e.g., also depending on defaults set in the main CMakeLists.txt)
if(NOT __processedUserDefaults)

  message(STATUS "Loading machine configuration for TACC Frontera supercomputer.\n"
    "This configuration has been tested (2021-07-19) using the following modules (with the default Intel stack loaded):\n"
    "$ module load boost hdf5\n"
    "and with Charm++ configured as follows:\n"
    "$ cmake -DNETWORK=mpi -DSMP=OFF -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort ..\n")

  set(CMAKE_CXX_COMPILER icpc CACHE STRING "")
  set(CMAKE_C_COMPILER icc CACHE STRING "")
  set(CMAKE_Fortran_COMPILER ifort CACHE STRING "")
  set(CMAKE_Fortran_FLAGS "-nofor-main" CACHE STRING "Default Fortran flags")

  # these flag(s) are currently only used when using openmp-simd optimizations
  # (to specify available/prefered instruction sets).
  # The Frontera User Guide recommends the following flag for their Cascade
  # Lake (CLX) Compute Nodes
  set(CONFIG_ARCH_FLAGS "-xCORE-AVX512")

  # Setting package paths (e.g., Grackle)

  # Mark done
  set(__processedUserDefaults ON)

else()

  if (USE_DOUBLE_PREC)
    string(APPEND CMAKE_Fortran_FLAGS " -real-size 64 -double-size 64")
  endif()

endif()
