cmake_minimum_required(VERSION 3.12)

# Machine config needs to be imported before project() to properly
# set compilers this way.
if (DEFINED Enzo-E_CONFIG)
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/config/${Enzo-E_CONFIG}.cmake)
    include(${CMAKE_CURRENT_SOURCE_DIR}/config/${Enzo-E_CONFIG}.cmake)
  else()
    message(FATAL_ERROR
      "Given machine configuration ${Enzo-E_CONFIG}.cmake not found in ${CMAKE_CURRENT_SOURCE_DIR}/config.")
  endif()
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
find_package(Grackle)
if (USE_GRACKLE)
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

# Define preprocessor definitions

## Precision defines
option(USE_DOUBLE_PREC "Use double precision. Turn off for single precision." ON)
if (USE_DOUBLE_PREC)
  add_compile_definitions(CONFIG_PRECISION_DOUBLE)
else()
  add_compile_definitions(CONFIG_PRECISION_SINGLE)
endif()
add_compile_definitions(SMALL_INTS)


option(new_output "Temporary setting for using new Output implementation" OFF)

set(node_size "64" CACHE STRING "Maximum number of procesess per shared-memory node (can be larger than needed)")
add_compile_definitions(CONFIG_NODE_SIZE=${node_size})
math(EXPR node_size_3 "${node_size} * 3")
add_compile_definitions(CONFIG_NODE_SIZE_3=${node_size_3})

option(trace "Print out detailed messages with the TRACE() series of statements" OFF)


option(verbose "Trace main phases" OFF)


option(trace_charm "Print out messages with the TRACE_CHARM() and TRACEPUP() series of statements" OFF)


option(debug "Whether to enable displaying messages with the DEBUG() series of \
  statements. Also writes messages to out.debug.<P> where P is the \
  (physical) process rank. Still requires the \"DEBUG\" group to be \
  enabled in Monitor (that is Monitor::is_active(\"DEBUG\") must be true for any output)" OFF)
option(debug_field "" OFF)
option(debug_field_face "" OFF)


option(check "Do extra run-time checking.  Useful for debugging, but can potentially slow calculations down" OFF)


option(debug_verbose "Print periodically all field values.  See src/Field/field_FieldBlock.cpp" OFF)

option(memory "Track dynamic memory statistics.  Can be useful, but can cause problems on some \
  systems that also override new [] () / delete [] ()" OFF)

option(balance "Enable charm++ dynamic load balancing" ON)
set(balancer "CommonLBs" CACHE STRING "Charm++ load balancer to use")


option(use_gprof "Compile with -pg to use gprof for performance profiling" OFF)


option(use_valgrind "Rn the test programs using valgrind to check for memory leaks" OFF)


option(use_performance "Use Cello Performance class for collecting performance \
  data (currently requires global reductions, and may not be fully \
  functional) (basic time data on root processor is still output)" ON)


option(use_projections "Compile the CHARM++ version for use with the Projections performance tool." OFF)


set(ip_charm "4" CACHE STRING "Number of processors to run parallel unit tests")


option(have_git "Is this is a Git repository" ON)


option (use_jemalloc "Use the jemalloc library for memory allocation" OFF)


# TODO(pgrete) this should be a package
option(use_papi "System has the PAPI performance API installed" OFF)

add_compile_definitions(NO_FREETYPE)
add_compile_definitions(CONFIG_USE_PERFORMANCE CONFIG_USE_MEMORY CONFIG_HAVE_VERSION_CONTROL)

find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(PNG REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)
# HDF5 Interface library
add_library(HDF5_C INTERFACE)
target_link_libraries(HDF5_C INTERFACE ${HDF5_C_LIBRARIES})
target_compile_definitions(HDF5_C INTERFACE ${HDF5_C_DEFINITIONS})
target_include_directories(HDF5_C INTERFACE ${HDF5_C_INCLUDE_DIRS})



set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src/Cello)
add_subdirectory(src/Enzo)
add_subdirectory(src/External)

# now generate the the config
configure_file(auto_config.def.in auto_config.def @ONLY)
message(STATUS
  "Cello/Enzo-E configured as follows:\n"
  "USE_GRACKLE = ${USE_GRACKLE}")
