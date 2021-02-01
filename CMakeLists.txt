
cmake_minimum_required(VERSION 3.12)

# Imports machine-specific configuration
#include(cmake/MachineCfg.cmake)

project(enzo-e LANGUAGES C CXX Fortran)

INCLUDE( FindPackageHandleStandardArgs )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Don't allow in-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR
      "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
      "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

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

#Get CHARM++ variables
# Checks an environment variable; note that the first check
# does not require the usual CMake $-sign.
IF( DEFINED ENV{CHARM_DIR} )
  SET( CHARM_DIR "$ENV{CHARM_DIR}" )
ELSE()
  message( STATUS "CHARM_DIR not defined")
ENDIF()


set(CHARMC_EXE "${CHARM_DIR}/bin/charmc")
include_directories(${CHARM_DIR}/include)

add_compile_definitions(CONFIG_PRECISION_DOUBLE)
add_compile_definitions(CONFIG_NODE_SIZE=4)
add_compile_definitions(NO_FREETYPE)


function(preprocess_charm_interfaces out_var)
  set(result)
  #Loop over each input file
  foreach(in_f ${ARGN})
    #Find the filename of the current file
    #file(RELATIVE_PATH out_f ${CMAKE_CURRENT_SOURCE_DIR} ${in_f})
    

	string(REGEX REPLACE "[.]CI" ".ci"  in_ci  ${in_f})
	configure_file(${in_f} ${in_ci} COPYONLY)

	#Create strings for the headers that charmc produces
    string(REGEX REPLACE "[.]ci" ".def.h"  out_def  ${in_ci})
    string(REGEX REPLACE "[.]ci" ".decl.h" out_decl ${in_ci})


    #set(out_f "${CMAKE_CURRENT_BINARY_DIR}/${out_def} " "${CMAKE_CURRENT_BINARY_DIR}/${out_decl} ")

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${out_def} ${CMAKE_CURRENT_BINARY_DIR}/${out_decl}
      COMMAND ${CHARMC_EXE}  ${in_ci} -I${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDS ${in_ci}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Creating charm++ headers ${out_f}"
      VERBATIM
      )
    #list(APPEND result ${out_f})
    list(APPEND result ${CMAKE_CURRENT_BINARY_DIR}/${out_def})
    list(APPEND result ${CMAKE_CURRENT_BINARY_DIR}/${out_decl})
  endforeach()
  #message(STATUS ${result} )
  set(${out_var} "${result}" PARENT_SCOPE)
endfunction()


add_subdirectory(src)

