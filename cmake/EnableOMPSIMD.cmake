# See LICENSE_ENZO file for license and copyright information

#.rst:
# EnableOMPSIMD
# -------
#
# Adds the openmp-simd flags to the C and C++ compilers (and any other flags
# necessary for simd optimizations).

if(__enableOMPSIMDModule)
  return()
endif()
set(__enableOMPSIMDModule YES)

# Function 'enableOMPSIMD' is used to add compiler flags for building with
# OpenMP's SIMD directives
# - For clang, gnu and intel compilers this does NOT enable other OpenMP
#   directives and should not link the openmp runtime library
function(enableOMPSIMD)

  # check assumption that CMAKE_C_COMPILER_ID and CMAKE_CXX_COMPILER_ID are
  # equal to each other...
  
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # theoretically includes Clang and AppleClang

    set(HOSTARCHFLAG "-march=native") # this is a hypothetical default option
    set(SPEED_FLAGS "-O3 -funroll-loops") # prioritize speed over size
    set(NEWFLAGS "-fopenmp-simd -ffast-math")

    message(FATAL_ERROR "enableOMPSIMD is untested for Clang.")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    # in case we ever want to support a version with OpenMP-SIMD, that does not
    # enable value-unsafe optimizations of floating-point operations, try the
    # following flags:
    # - "-fopenmp-simd -funroll-loops -fno-trapping-math -fno-signed-zeros"
    # - if we also include -fno-math-errno, it will cause 1D MHD shock tube
    #   problems to have slightly different L1 error norms for the VL+CT
    #   solver, depending on the position in the transverse direction
    #   (see MHD_shock_tube_test)

    set(HOSTARCHFLAG "-march=native") # this is a hypothetical default option
    set(SPEED_FLAGS "-O3 -funroll-loops") # prioritize speed over size
    set(NEWFLAGS "-fopenmp-simd -ffast-math")

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    # this matches both the legacy icc and icpc compiler...
    # For CMake Ver >= 20, CMAKE_CXX_COMPILER_ID is "IntelLLVM" for the latter
    # (but this will continue to match both)

    # technically -qopenmp-simd is already enabled by default when compiling
    # with -O2 or higher
    set(HOSTARCHFLAG "-xHost") # this is a hypothetical default option
    set(SPEED_FLAGS "-O3") # prioritize speed over size
    set(NEWFLAGS "-qopenmp-simd")

  else()
    message(FATAL_ERROR
      "OpenMP-SIMD handling is not yet implemented for the "
      "${CMAKE_CXX_COMPILER_ID} compiler."
      )
  endif()

  if (NOT DEFINED CONFIG_ARCH_FLAGS)
    message(FATAL_ERROR
      "The CONFIG_ARCH_FLAGS variable is not defined.\n"
      "This variable is strongly-recommended while compiling with OpenMP-SIMD, in order to inform the compiler which vector instructions are available/prefered.\n"
      "This variable should be defined in the machine config file.\n"
      "A reasonable default value for this variable, when using the ${CMAKE_CXX_COMPILER_ID} compiler, might be \"${HOSTARCHFLAG}\" (this tells the compiler to target the CPU architecture of the machine that is performing the compilation)."
      )
  endif()

  # now update the flags!
  # - we NEED to update CMAKE_<LANG>_FLAGS_<CONFIG>
  # - if we update CMAKE_<LANG>_FLAGS the optimization level flags in
  #   CMAKE_<LANG>_FLAGS_<CONFIG> get precedence. This is specifically an
  #   issue for RelWithDebInfo builds
  # - we could theoretically use add_compile_options(), but that would also
  #   pass these flags to the FORTRAN compiler
  foreach(BTYPE IN ITEMS "DEBUG" "RELEASE" "MINSIZEREL" "RELWITHDEBINFO")
    foreach(LANG IN ITEMS "C" "CXX")

      # make a copy of global variable to be modified
      set(localCopy ${CMAKE_${LANG}_FLAGS_${BTYPE}})

      # append the new flags to the local copy
      string(APPEND localCopy " ${NEWFLAGS} ${CONFIG_ARCH_FLAGS}")

      if ((BTYPE STREQUAL "RELEASE") OR
          (BTYPE STREQUAL "RELWITHDEBINFO"))
        # append strings that prioritize speed (at expense of size)
        string(APPEND localCopy " ${SPEED_FLAGS}")
      endif()

      # overwrite the global variable
      set(CMAKE_${LANG}_FLAGS_${BTYPE} ${localCopy} PARENT_SCOPE)

    endforeach(LANG)
  endforeach(BTYPE)
endfunction(enableOMPSIMD)

# REQUIRE TARGET_ARCH_FLAGS for gcc compiler...
