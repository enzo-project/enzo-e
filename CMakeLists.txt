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
# TODO set proper version before merge
project(Enzo-E VERSION 0.9.0 LANGUAGES C CXX Fortran)

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

find_package(Charm REQUIRED)
# Link executables with the charmc wrapper
STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>" "${CHARM_COMPILER}"
       CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")

# Need to process Grackle here as both Cello and Enzo-E depend on it
option(USE_GRACKLE "Use Grackle Chemistry" ON)
if (USE_GRACKLE)
  find_package(Grackle)
  if (Grackle_FOUND)
    # Setting global compile def as performance and simulation Cello libs need the var
    # We should discuss the use of a central config.hpp (similar to the existing
    # auto_config.def but also centrally included beyond informational output).
    add_compile_definitions(CONFIG_USE_GRACKLE)

    # Also specifically setting the defines used in the *.ci files of Cello and Enzo.
    # Alternatively, figure out a way to extract those from the `COMPILE_DEFINITIONS`
    # property and directly process in `cmake/charm.cmake`.
    set(CHARM_PREPROC_DEFS ${CHARM_PREPROC_DEFS} "-DCONFIG_USE_GRACKLE ")
  else()
    message(FATAL_ERROR
      "Configured to use Grackle but Grackle library not found.\n"
      "Either disable grackle (e.g., `-DUSE_GRACKLE=OFF`) or provide path "
      "(e.g., `-DGrackle_ROOT=/PATH/TO/GRACKLE/INSTALL/DIRECTORY`).")
  endif()
endif()

find_package(Boost REQUIRED COMPONENTS filesystem)
# Fix linking error to
# /opt/apps/gcc/8.3.0/bin/ld: main.cpp:(.text+0x679e): undefined reference to `boost::system::generic_category()'
# /opt/apps/gcc/8.3.0/bin/ld: main.cpp:(.text+0x67aa): undefined reference to `boost::system::system_category()'
# due to deprecated code (`filesystem` depends on `system`).
add_compile_definitions(BOOST_SYSTEM_NO_DEPRECATED)

find_package(PNG REQUIRED)
add_compile_definitions(NO_FREETYPE)

find_package(HDF5 REQUIRED COMPONENTS C)
# HDF5 Interface library
add_library(HDF5_C INTERFACE)
target_link_libraries(HDF5_C INTERFACE ${HDF5_C_LIBRARIES})
target_compile_definitions(HDF5_C INTERFACE ${HDF5_C_DEFINITIONS})
target_include_directories(HDF5_C INTERFACE ${HDF5_C_INCLUDE_DIRS})


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


option (use_jemalloc "Use the jemalloc library for memory allocation" OFF)
if (use_jemalloc)
  find_package(jemalloc)
  if (jemalloc_FOUND)
    add_compile_definitions(CONFIG_USE_JEMALLOC)
  else()
    message(FATAL_ERROR
      "Requested to use the jemalloc library for memory allocation but jemalloc was not found. "
      "Try setting specific path via `-Djemalloc_ROOT=/PATH/TO/jemalloc/INSTALL` "
      " or disable jemalloc via `-Duse_jemalloc=OFF` (default)."
      )
  endif()
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


option(use_papi "Use the PAPI performance API" OFF)
if (use_papi)
  find_package(PAPI)
  if (PAPI_FOUND)
    add_compile_definitions(CONFIG_USE_PAPI PAPI3)
  else()
    message(FATAL_ERROR
      "Requested to use PAPI performance API but PAPI was not found. "
      "Try setting specific path via `-DPAPI_ROOT=/PATH/TO/PAPI/INSTALL` "
      " or disable PAPI via `-Duse_papi=OFF` (default)."
      )
  endif()
endif()

option(USE_SIMD "Use OpenMP SIMD directives. Also enables value-unsafe floating-point optimizations which may already be enable on some compilers (like icc)." OFF)
include("cmake/EnableOMPSIMD.cmake")
if (USE_SIMD)
  enableOMPSIMD()
endif()

# Include machine file second time (if used initially) to set additional options that may depend
# on global (default) options, such as USE_DOUBLE_PREC
if (__processedUserDefaults)
  include(${Enzo-E_CONFIG_PATH})
endif()

# In principle a more fine grained control (i.e., target specific include directories
# rather than this global one would be preferred, but, e.g., `pngwriter.h` is curently
# included in many targets/libraries (without being linked) so we'll use the global for
# convenience for now.
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/External)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src/Cello)
add_subdirectory(src/Enzo)
add_subdirectory(src/External)

include(CTest)

# including the CTest module defines the BUILD_TESTING option. Its initial value is ON
option(USE_YT_BASED_TESTS "Includes tests that require yt" ON)
if (BUILD_TESTING)
  find_package(Python3)
  if (Python3_FOUND)
    message(STATUS
      "   IMPORTANT: Make sure this is the Python version you want "
      "to use. If not, e.g., because you want to use a conda or virtualenv based "
      "version, try setting `-DPython3_FIND_VIRTUALENV=ONLY` or "
      "`-DPython3_EXECUTABLE=/path/to/desired/python`."
      )
  else()
    message(FATAL_ERROR
     "You asked to build the tests, which require a python3 installation, but no python "
     "installation was found. You can either install a python distribution (recommended) or "
     "disable all tests by setting `-DBUILD_TESTING=OFF`."
     )
  endif()

  if (USE_YT_BASED_TESTS)
    # Ensure all required packages are present
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/PythonModuleCheck.cmake)
    required_python_modules_found("yt;libconf;h5py")
  endif()

  add_subdirectory(test)
endif()


# extract compile defs from from Cello to populate config
get_directory_property( CELLO_CPPDEFINES DIRECTORY src/Cello COMPILE_DEFINITIONS )
# now generate the the config
configure_file(auto_config.def.in auto_config.def ESCAPE_QUOTES @ONLY)
