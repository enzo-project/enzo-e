
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

LINK_DIRECTORIES(${CHARM_DIR}/lib)
set(CHARM_LIBRARIES "moduleCommonLBs" "ckmain" "ck" "memory-default" "threads-default" "conv-machine" "conv-core" "tmgr" "conv-util" "conv-partition" "conv-ldb" "memory-default" "threads-default" "ldb-rand" "ckqt" "moduleCommonLBs" "moduleNDMeshStreamer" "modulecompletion" "hwloc_embedded" "mpi-mainmodule")
set(ENZO_E_LIBRARIES  "dl" "pthread" "z" "hdf5" "png" "gfortran")

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost REQUIRED COMPONENTS system filesystem) 

#Cello variables

add_compile_definitions(CONFIG_PRECISION_DOUBLE)
add_compile_definitions(CONFIG_NODE_SIZE=4)
add_compile_definitions(CONFIG_NODE_SIZE_3=12)
add_compile_definitions(NO_FREETYPE)
add_compile_definitions(SMALL_INTS)

#set(CHARMC_EXE "${CHARM_DIR}/bin/charmc")

#define CELLO_ARCH "linux_gnu"
#define CELLO_PREC "double"
#define CELLO_CC "gcc "
#define CELLO_CFLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic   "
set( CMAKE_C_FLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic   ")
#define CELLO_CPPDEFINES "CONFIG_PRECISION_DOUBLE SMALL_INTS {'CONFIG_NODE_SIZE': 64} {'CONFIG_NODE_SIZE_3': 192} NO_FREETYPE CONFIG_USE_PERFORMANCE CONFIG_USE_MEMORY CONFIG_HAVE_VERSION_CONTROL"
#define CELLO_CPPPATH "#/include /usr/include /usr/include/boost /lib/x86_64-linux-gnu/include"
#define CELLO_CXX "/home/forrest/MSU/research/code/charm//bin/charmc -language charm++  "
#define CELLO_CXXFLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic     -balancer CommonLBs"
set( CMAKE_CXX_FLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic" )
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffast-math -funroll-loops -fPIC -pedantic -ffixed-line-length-132    -fdefault-real-8 -fdefault-double-8")
#define CELLO_FORTRAN "gfortran"
#define CELLO_FORTRANLIBS "gfortran"
#define CELLO_FORTRANPATH "#/include"
#define CELLO_LIBPATH "#/lib /usr/lib /usr/lib/x86_64-linux-gnu /lib/x86_64-linux-gnu/lib /usr/lib/x86_64-linux-gnu"
#define CELLO_LINKFLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic     -module CommonLBs"
set( CMAKE_EXE_LINKER_FLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic -module CommonLBs" )
#define CELLO_HOST "forrest-razer"
#define CELLO_DIR "/home/forrest/MSU/research/code/enzo-e"
#define CELLO_DATE "2021-02-02"
#define CELLO_TIME "00:23:53"
#define CHARM_VERSION 61002
#define CHARM_PATH "/home/forrest/MSU/research/code/charm/"
#define CELLO_CHANGESET "f04dd0621a4e63e2c39d8390f4b93f2c57e6a3f1"
#define CHARM_BUILD "multicore-linux-x86_64"

#configure_file(${CMAKE_SOURCE_DIR}/auto_config.def.in ${CMAKE_BINARY_DIR}/auto_config.def)

message( STATUS ${CHARMC_EXE})
#SET(CMAKE_C_LINK_EXECUTABLE ${CHARMC_EXE})
#SET(CMAKE_CXX_LINK_EXECUTABLE ${CHARMC_EXE})
SET(CMAKE_LINKER ${CHARMC_EXE})
set(CMAKE_CXX_LINK_EXECUTABLE "<CMAKE_LINKER> ${CMAKE_EXE_LINKER_FLAGS} <OBJECTS> <LINK_LIBRARIES> -o <TARGET>")

function(preprocess_charm_interfaces out_var)
  set(result)
  #Loop over each input file
  foreach(in_f ${ARGN})
    #Find the filename of the current file

	string(REGEX REPLACE "[.]CI" ".ci"  in_ci  ${in_f})
	configure_file(${in_f} ${in_ci} COPYONLY)

	#Create strings for the headers that charmc produces
    string(REGEX REPLACE "[.]ci" ".def.h"  out_def  ${in_ci})
    string(REGEX REPLACE "[.]ci" ".decl.h" out_decl ${in_ci})
    #Remove directory name (HACK?)
    get_filename_component(out_def  ${out_def}  NAME)
    get_filename_component(out_decl ${out_decl} NAME)


    #set(out_f "${CMAKE_CURRENT_BINARY_DIR}/${out_def} " "${CMAKE_CURRENT_BINARY_DIR}/${out_decl} ")

    message(STATUS "${in_f} ${out_def} ${out_decl}")

    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${out_def} ${CMAKE_CURRENT_BINARY_DIR}/${out_decl}
      COMMAND ${CHARMC_EXE}  ${in_ci} -I${CMAKE_SOURCE_DIR}/src/Cello -I${CMAKE_SOURCE_DIR}/src/Enzo
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${in_f}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Creating charm++ headers ${out_def} ${out_decl}"
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

