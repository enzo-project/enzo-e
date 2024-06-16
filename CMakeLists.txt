cmake_minimum_required(VERSION 3.16)

# Machine config needs to be imported before project() to properly
# set compilers this way.
if (DEFINED Enzo-E_CONFIG)
  if (EXISTS $ENV{HOME}/.enzo-e/${Enzo-E_CONFIG}.cmake)
    set(Enzo-E_CONFIG_PATH $ENV{HOME}/.enzo-e/${Enzo-E_CONFIG}.cmake)
  elseif (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/config/${Enzo-E_CONFIG}.cmake)
    set(Enzo-E_CONFIG_PATH ${CMAKE_CURRENT_SOURCE_DIR}/config/${Enzo-E_CONFIG}.cmake)
  else()
    message(FATAL_ERROR
      "Given machine configuration ${Enzo-E_CONFIG}.cmake neither found in $ENV{HOME}/.enzo-e "
      "nor in ${CMAKE_CURRENT_SOURCE_DIR}/config.")
  endif()
  message(STATUS "Using machine configuration file from ${Enzo-E_CONFIG_PATH}")
  include(${Enzo-E_CONFIG_PATH})
endif()

#Define project and languages
project(Enzo-E VERSION 1.0.0 LANGUAGES C CXX Fortran)

# We need C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Don't allow in-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR
      "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
      "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Ensure the custom modules to find Charm++ and Grackle are included
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# confirm that std::filesystem is supported by the compiler and create an
# interface target called StdFilesystem::StdFilesystem that encapsulates
# usage requirements (if there are any)
include("CreateStdFilesystemTarget")
create_StdFilesystem_target()

#----------------------------------------------------------------------------------------
# External libs

# If the user doesn't specify a build type, prefer Release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include("dependencies.cmake")

set(Cello_TARGET_LINK_OPTIONS "")

#----------------------------------------------------------------------------------------
# Define preprocessor definitions/user configuration options

option(USE_DOUBLE_PREC "Use double precision. Turn off for single precision." ON)
if (USE_DOUBLE_PREC)
  add_compile_definitions(CONFIG_PRECISION_DOUBLE)
  # PREC_STRING is used during test setup
  set(PREC_STRING "double")
else()
  add_compile_definitions(CONFIG_PRECISION_SINGLE)
  # PREC_STRING is used during test setup
  set(PREC_STRING "single")
endif()
add_compile_definitions(SMALL_INTS)

add_compile_definitions(CELLO_VERSION="${CMAKE_PROJECT_VERSION}")

# Whether to bypass passing MsgRefine directly to Block constructor,
# or request it from a separate entry method to bypass a Charm++
# memory leak. This should only be set to 0 after (and if) the bug is
# addressed in Charm++, or when explicitly testing a Charm++ build for
# this bug.
option(bypass_charm_mem_leak "Temporary setting to bypass Charm++ bug" ON)
if (bypass_charm_mem_leak)
  add_compile_definitions(BYPASS_CHARM_MEM_LEAK)
  set(CHARM_PREPROC_DEFS ${CHARM_PREPROC_DEFS} "-DBYPASS_CHARM_MEM_LEAK ")
endif()

set(node_size "64" CACHE STRING "Maximum number of procesess per shared-memory node (can be larger than needed)")
add_compile_definitions(CONFIG_NODE_SIZE=${node_size})
math(EXPR node_size_3 "${node_size} * 3")
add_compile_definitions(CONFIG_NODE_SIZE_3=${node_size_3})

option(trace "Print out detailed messages with the TRACE() series of statements" OFF)
if (trace)
  add_compile_definitions(CELLO_TRACE)
endif()

option(verbose "Trace main phases" OFF)
if (verbose)
  add_compile_definitions(CELLO_VERBOSE)
endif()


option(trace_charm "Print out messages with the TRACE_CHARM() and TRACEPUP() series of statements" OFF)
if (trace_charm)
  add_compile_definitions(CELLO_TRACE_CHARM)
endif()


option(debug "Whether to enable displaying messages with the DEBUG() series of \
  statements. Also writes messages to out.debug.<P> where P is the \
  (physical) process rank. Still requires the \"DEBUG\" group to be \
  enabled in Monitor (that is Monitor::is_active(\"DEBUG\") must be true for any output)" OFF)
option(debug_field "" OFF)
option(debug_field_face "" OFF)
if (debug)
  add_compile_definitions(CELLO_DEBUG)
endif()
if (debug_field)
  add_compile_definitions(DEBUG_FIELD)
endif()
if (debug_field_face)
  add_compile_definitions(DEBUG_FIELD_FACE)
endif()


option(check "Do extra run-time checking.  Useful for debugging, but can potentially slow calculations down" OFF)
if (check)
  add_compile_definitions(CELLO_CHECK)
endif()


option(debug_verbose "Print periodically all field values.  See src/Field/field_FieldBlock.cpp" OFF)
if (debug_verbose)
  add_compile_definitions(CELLO_DEBUG_VERBOSE)
endif()

option(memory "Track dynamic memory statistics.  Can be useful, but can cause problems on some \
  systems that also override new [] () / delete [] ()" OFF)
if (memory)
  add_compile_definitions(CONFIG_USE_MEMORY)
endif()

option(balance "Enable charm++ dynamic load balancing" ON)
set(balancer_included "CommonLBs" CACHE STRING "Charm++ load balancer to include")
set(balancer_default "TreeLB" CACHE STRING "Charm++ load balancer to use by default")
#TODO we should figure our reasonale defaults and/or provide instructions in the docs
if (balance)
  foreach(BALANCER IN LISTS balancer_included)
    string(APPEND Cello_TARGET_LINK_OPTIONS " -module ${BALANCER}")
  endforeach()
  foreach(BALANCER IN LISTS balancer_default)
    string(APPEND Cello_TARGET_LINK_OPTIONS " -balancer ${BALANCER}")
  endforeach()
endif()


option(use_gprof "Compile with -pg to use gprof for performance profiling" OFF)
if (use_gprof)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -pg")
endif()


option(use_performance "Use Cello Performance class for collecting performance \
  data (currently requires global reductions, and may not be fully \
  functional) (basic time data on root processor is still output)" ON)
if (use_performance)
  add_compile_definitions(CONFIG_USE_PERFORMANCE)
endif()


option(use_projections "Compile the CHARM++ version for use with the Projections performance tool." OFF)
if (use_projections)
  add_compile_definitions(CONFIG_USE_PROJECTIONS)
  string(APPEND Cello_TARGET_LINK_OPTIONS " -tracemode projections")
endif()


set(ip_charm "4" CACHE STRING "Number of processors to run parallel unit tests (legacy name)")
set(PARALLEL_LAUNCHER ${CHARM_RUN} CACHE STRING "Launcher to use for parallel tests")
set(PARALLEL_LAUNCHER_NPROC_ARG "+p" CACHE STRING "Argument to set number of processing elements for parallel launcher")
set(PARALLEL_LAUNCHER_NPROC ${ip_charm} CACHE STRING "Number of processors to run parallel unit tests")
if ("${PARALLEL_LAUNCHER_NPROC_ARG}" MATCHES "[ \t\r\n]")
  message(FATAL_ERROR
    "Whitespace should NOT be present in PARALLEL_LAUNCHER_NPROC_ARG.\n"
    "If you want to use PARALLEL_LAUNCHER_NPROC_ARG to pass multiple "
    "arguments to the launcher, use a semicolon to delimit each argument.\n"
    "For example, instead of setting PARALLEL_LAUNCHER_NPROC_ARG to "
    "\"++local +p\", set it equal to \"++local;+p\".")
endif()


option(have_git "Is this a Git repository" ON)
if (have_git)
  add_compile_definitions(CONFIG_HAVE_VERSION_CONTROL)
  # Extract current changeset.
  # Note, this will only be called during initial configure so changes
  # between configure and build may result in an outdated hash.
  execute_process(
      COMMAND
        git rev-parse --short HEAD
          RESULT_VARIABLE CHANGESET_RESULT
          OUTPUT_VARIABLE CELLO_CHANGESET
          OUTPUT_STRIP_TRAILING_WHITESPACE
        )
endif()

option (smp "Use Charm++ in SMP mode." OFF)
if (smp)
  if (CHARM_SMP)
    add_compile_definitions(CONFIG_SMP_MODE)
  else()
    message(FATAL_ERROR
      "Requested to use SMP in Cello/Enzo-E but could not find SMP support in Charm++. "
      "Either recompile Charm++ with SMP support or set `-Dsmp=OFF` (default) in Cello/Enzo-E."
      )
  endif()
endif()


# define recipies for building external dependencies before we introduce
# compiler flags specific to Enzo-E (and Cello)
add_subdirectory(src/External)

# fetch known sets of compiler flags
# ----------------------------------
# right now, we apply these flags with add_compile_options
# -> this applies the flags to when compiling all targets defined after this
#    point in the cmake-file (aside: add_compile_options behaves a little more
#    sensibly than other commands that apply global properties)
# -> to selectively apply flags based on language (and build-type), we
#    currently make use of generator expressions.
#
# NOTE: we should transition to defining a series of INTERFACE libraries
# (e.g. CommonOptionsC, CommonOptionsCXX, CommonOptionsFortran) and use
# those to propogate our compiler-flags to our various targets. This is
# desirable for 2 reasons:
#   1. it is a little more explicit and idiomatic
#   2. language-specific generator-expressions are not supported on all
#      build-generators (this isn't really a problem per-se, but it would be
#      nice to be as flexible as possible)

# handle flags within ENZOE_<LANG>_FLIST for each language used be Enzo-E
include("GenericOptionCommand")
get_property(_ENZOE_LANG_LIST GLOBAL PROPERTY ENABLED_LANGUAGES)

foreach(lang IN LISTS _ENZOE_LANG_LIST)
  # first, initialize ENZOE_<LANG>_FLIST if it hasn't already been initialized
  generic_option(
    ENZOE_${lang}_FLIST STRING "\
specifies flags that are passed to the compiler (in all configurations) for \
compiling the ${lang} files that are used in Cello and Enzo-E. There are 3 \
main differences from CMAKE_${lang}_FLAGS: \
(i) these flags won't be passed to any external dependencies that are \
compiled as part of this build \
(ii) these flags will NOT be passed to the linker, and \
(ii) these flags must be specified as a semi-colon delimited list of values."
    # if there are any, the default values would be set in machine file
    DEFAULT_VALUE "${ENZOE_${lang}_FLIST_INIT}"
    FORBID_EXISTING_NONCACHE_VAR # <- discourage bad habits!
  )

  # now setup cmake to make use of any compiler flags in the list
  if(NOT ("${ENZOE_${LANG}_FLIST}" STREQUAL ""))
    message("$<$<COMPILE_LANGUAGE:${LANG}>:${ENZOE_${LANG}_FLIST}>")
  endif()
endforeach()

# introduce floating-point related compiler options
option(OPTIMIZE_FP "Enables value-unsafe floating-point optimizations (this may already be enabled on some compilers like icc)." OFF)
option(USE_SIMD "Use OpenMP SIMD directives. This may already be enabled on some compilers (like icc)." OFF)
include("FPOptimizationCompileOptions")
if (OPTIMIZE_FP)
  get_fp_optimization_compile_options(${USE_SIMD} _ENZOE_optimization_options)
  add_compile_options("${_ENZOE_optimization_options}")
elseif (USE_SIMD)
  message(FATAL_ERROR "Can only use `USE_SIMD=ON` when `OPTIMIZE_FP=ON`.")
endif()

# introduce required Fortran Compiler Options
include("RequiredFortranCompileOptions")
get_required_fortran_options(_ENZOE_REQ_FORTRAN_OPTS)
add_compile_options("${_ENZOE_REQ_FORTRAN_OPTS}")

# Include machine file second time (if used initially) to set additional options that may depend
# on global (default) options, such as USE_DOUBLE_PREC
if (__processedUserDefaults)
  include(${Enzo-E_CONFIG_PATH})
endif()

# we need to pull CTest sooner or later anyways. We do it before adding the
# source directories because it introduces the BUILD_TESTING option which is
# set to ON if it was not previously initialized (e.g. from the command line)
# -> this option is used within source directories to determine whether the
#    binaries for unit-tests should be declared
include(CTest)

# In principle a more fine grained control (i.e., target specific include directories
# rather than this global one would be preferred, but, e.g., `pngwriter.h` is curently
# included in many targets/libraries (without being linked) so we'll use the global for
# convenience for now.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/External)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(USE_PRECOMPILED_HEADERS
  "Use precompiled headers to speed up compilation" ON)

add_subdirectory(src/Cello)
add_subdirectory(src/Enzo)


if (BUILD_TESTING)
  add_subdirectory(test)
endif()

# extract compile defs from from Cello to populate config
get_directory_property( CELLO_CPPDEFINES DIRECTORY src/Cello COMPILE_DEFINITIONS )
# now generate the the config
configure_file(auto_config.def.in auto_config.def ESCAPE_QUOTES @ONLY)
