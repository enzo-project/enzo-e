# See LICENSE_ENZO file for license and copyright information

#.rst:
# EnableFPOptimizations
# ---------------------
#
# Defines a function that enables value-unsafe floating point optimization
# flags for the C and C++ compilers.
# 
# This can also add the openmp-simd flags to the C and C++ compilers (and any
# other flags necessary for simd optimizations).

if(__enableFPOptimizations)
  return()
endif()
set(__enableFPOptimizations YES)

# Function 'enableFPOptimizations' is used to add C and C++ flags to the
# various build types that compile the program with value unsafe floating point
# optimizations.
#
# ARGUMENTS
# ---------
# USE_SIMD
#   When this a passed a TRUE value, this function also adds C and C++ compiler
#   flags to enable OpenMP's SIMD directives
#   - For clang, gnu and intel compilers this does NOT enable other OpenMP
#     directives and should not link the openmp runtime library
#
# NOTES
# -----
# This function assumes that the global CONFIG_ARCH_FLAGS variable has been
# defined by the machine configuration file (this flag should provide the
# compiler about the architecture of the CPU where the code will be executed,
# which lets the compiler choose optimal instructions). If this variable has
# not be defined, the function will produce an error and provide the user with
# a sensible default that they may want to use.
#
# This function may have less impact on the performance of programs compiled by
# the Intel compiler compared to code compiled by other compilers. This is
# because the Intel compilers enable some of these options by default (this
# strongly to the Intel Compiler's reputation for producing faster code)
function(enableFPOptimizations USE_SIMD)

  # ToDo: check assumption that CMAKE_C_COMPILER_ID and CMAKE_CXX_COMPILER_ID
  # are equal to each other...

  # First, determine the values for the following variables based on the
  # compiler type:
  # - DFLT_HOSTARCHFLAG:
  #     * this stores the flag telling the compiler that it can assume that the
  #       CPU architecture of the machine currently performing the compilation
  #       is identical to that of the machine where the compiled program is run
  #     * this variable is only used to describe a potential default value that
  #       the user could use in the error message raised when the
  #       CONFIG_ARCH_FLAGS global variable was not defined in the machine
  #       configuration file.
  # - SPEED_FLAGS:
  #     * Specifies flags that prioritize code speed over code size.
  #     * Only the "RELEASE" and "RELWITHDEBINFO" build types will make use of
  #       these flags
  # - FP_FLAGS:
  #     * Specifies flags enabling value-unsafe floating point optimizations.
  # - OMPSIMD_FLAGS:
  #     * Specifies flags enabling OpenMP's SIMD directives

  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    # theoretically includes Clang and AppleClang

    set(DFLT_HOSTARCHFLAG "-march=native")
    set(SPEED_FLAGS "-O3 -funroll-loops")
    set(FP_FLAGS "-fopenmp-simd")
    set(OMPSIMD_FLAGS "-fopenmp-simd")

    message(FATAL_ERROR "enableOFPOptimizations is untested for Clang.")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

    # in case we ever want to support a version with OpenMP-SIMD, that does not
    # enable value-unsafe optimizations of floating-point operations, try the
    # following flags:
    # - "-fopenmp-simd -funroll-loops -fno-trapping-math -fno-signed-zeros"
    # - if we also include -fno-math-errno, it will cause 1D MHD shock tube
    #   problems to have slightly different L1 error norms for the VL+CT
    #   solver, depending on the position in the transverse direction
    #   (see MHD_shock_tube_test)

    set(DFLT_HOSTARCHFLAG "-march=native")
    set(SPEED_FLAGS "-O3 -funroll-loops")
    set(FP_FLAGS "-ffast-math")
    set(OMPSIMD_FLAGS "-fopenmp-simd")

  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    # this matches both the legacy icc and icpc compiler...
    # For CMake Ver >= 20, CMAKE_CXX_COMPILER_ID is "IntelLLVM" for the latter
    # (but this will continue to match both)

    # technically -qopenmp-simd is already enabled by default when compiling
    # with -O2 or higher
    set(DFLT_HOSTARCHFLAG "-xHost")
    set(SPEED_FLAGS "-O3")
    set(FP_FLAGS "")
    set(OMPSIMD_FLAGS "-qopenmp-simd")

  else()
    message(FATAL_ERROR
      "OpenMP-SIMD handling is not yet implemented for the "
      "${CMAKE_CXX_COMPILER_ID} compiler."
      )
  endif()

  # Second, check if the CONFIG_ARCH_FLAGS variable is enabled
  if (NOT DEFINED CONFIG_ARCH_FLAGS)
    message(FATAL_ERROR
      "The CONFIG_ARCH_FLAGS variable is not defined.\n"
      "This variable is strongly-recommended while compiling with OpenMP-SIMD, in order to inform the compiler which vector instructions are available/prefered.\n"
      "This variable should be defined in the machine config file.\n"
      "A reasonable default value for this variable, when using the ${CMAKE_CXX_COMPILER_ID} compiler, might be \"${HOSTARCHFLAG}\" (this tells the compiler to target the CPU architecture of the machine that is performing the compilation)."
      )
  endif()

  # Finally, update the flags
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
      string(APPEND localCopy " ${FP_FLAGS} ${CONFIG_ARCH_FLAGS}")

      if(USE_SIMD)
        string(APPEND localCopy " ${OMPSIMD_FLAGS}")
      endif()

      if ((BTYPE STREQUAL "RELEASE") OR
          (BTYPE STREQUAL "RELWITHDEBINFO"))
        # append strings that prioritize speed (at expense of size)
        string(APPEND localCopy " ${SPEED_FLAGS}")
      endif()

      # overwrite the global variable
      set(CMAKE_${LANG}_FLAGS_${BTYPE} ${localCopy} PARENT_SCOPE)

    endforeach(LANG)
  endforeach(BTYPE)
endfunction(enableFPOptimizations)

# REQUIRE TARGET_ARCH_FLAGS for gcc compiler...
