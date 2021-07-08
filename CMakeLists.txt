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

#Define preprocessor definitions
add_compile_definitions(CONFIG_PRECISION_DOUBLE)
add_compile_definitions(CONFIG_NODE_SIZE=4)
add_compile_definitions(CONFIG_NODE_SIZE_3=12)
add_compile_definitions(NO_FREETYPE)
add_compile_definitions(SMALL_INTS)
add_compile_definitions(CONFIG_USE_PERFORMANCE CONFIG_USE_MEMORY CONFIG_HAVE_VERSION_CONTROL)

find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(PNG REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(src/Cello)
add_subdirectory(src/Enzo)
add_subdirectory(src/External)

message(STATUS
  "Cello/Enzo-E configured as follows:\n"
  "USE_GRACKLE = ${USE_GRACKLE}")
