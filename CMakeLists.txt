cmake_minimum_required(VERSION 3.12)

#Define project and languages
# TODO set proper version before merge
project(enzo-e VERSION 0.9.0 LANGUAGES C CXX Fortran)

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

