# See LICENSE_ENZO file for license and copyright information

#.rst:
# FPOptimizationCompileOptions
# ----------------------------
#
# Defines a function that enables value-unsafe floating point optimization
# flags for the C and C++ compilers.
#
# This can also add the openmp-simd flags to the C and C++ compilers (and any
# other flags necessary for simd optimizations).

if(__FPOptimizationCompileOptions)
  return()
endif()
set(__FPOptimizationCompileOptions YES)

# Function 'get_fp_optimization_compile_options' is used to retrieve C and C++
# flags, that compile the program with value unsafe floating point
# optimizations. The flags use generator-expressions to conditionally
# enable/disable choices based on the precise built-type
#
# ARGUMENTS
# ---------
# USE_SIMD
#   When this is passed a TRUE value, this function also adds C and C++
#   compiler flags to enable OpenMP's SIMD directives
#   - For clang, gnu and intel compilers this does NOT enable other OpenMP
#     directives and should not link the openmp runtime library
# outVar
#   The name of the variable where the list of compiler-options are written to
#   by this function.
#
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
function(get_fp_optimization_compile_options USE_SIMD outVar)

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

  if(CMAKE_CXX_COMPILER_ID MATCHES "^AppleClang|Clang$")

    set(DFLT_HOSTARCHFLAG "-march=native")
    set(SPEED_FLAGS "-O3;-funroll-loops")
    set(FP_FLAGS "-ffast-math")

    if (USE_SIMD AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
      message(FATAL_ERROR "AppleClang does not support OpenMP-SIMD")
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set(OMPSIMD_FLAGS "-fopenmp-simd")
      message(FATAL_ERROR
        "get_fp_optimization_compile_options is untested for Clang.")
    endif()

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
    set(SPEED_FLAGS "-O3;-funroll-loops")
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
      "get_fp_optimization_compile_options does not support the "
      "${CMAKE_CXX_COMPILER_ID} compiler yet."
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

  # Finally, construct the list of flags
  set(optionList ${CONFIG_ARG_FLAGS} ${FP_FLAGS})

  if(USE_SIMD)
    list(APPEND optionList ${OMPSIMD_FLAGS})
  endif()

  if (NOT (SPEED_FLAGS STREQUAL ""))
    # ONLY include SPEED_FLAGS if building RELEASE or RELWITHDEBINFO configs
    # -> these flags prioritize speed at expense of size and debugability
    # -> starting in cmake version 3.19, could be more concise
    list(APPEND optionList
      "$<$<CONFIG:RELEASE>:${SPEED_FLAGS}>"
      "$<$<CONFIG:RELWITHDEBINFO>:${SPEED_FLAGS}>")
  endif()

  set("${outVar}" "$<$<COMPILE_LANGUAGE:C,CXX>:${optionList}>" PARENT_SCOPE)

endfunction()
