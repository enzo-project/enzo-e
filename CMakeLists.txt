cmake_minimum_required(VERSION 3.12)

project(enzo-e VERSION 0.9.0 LANGUAGES C CXX Fortran)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

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

