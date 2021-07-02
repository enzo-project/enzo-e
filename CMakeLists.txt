
cmake_minimum_required(VERSION 3.12)

#Import cmake package handling
INCLUDE( FindPackageHandleStandardArgs )


#Define project and languages
project(enzo-e LANGUAGES C CXX Fortran)


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
IF( DEFINED ENV{CHARM_HOME} )
  SET( CHARM_HOME "$ENV{CHARM_HOME}" )
ELSE()
  message( STATUS "CHARM_HOME not defined")
ENDIF()
set(CHARMC_EXE "${CHARM_HOME}/bin/charmc")
include_directories(${CHARM_HOME}/include)

#LINK_DIRECTORIES(${CHARM_HOME}/lib)
#set(CHARM_LIBRARIES "moduleCommonLBs" "ckmain" "ck" "memory-default" "threads-default" "conv-machine" "conv-core" "tmgr" "conv-util" "conv-partition" "conv-ldb" "memory-default" "threads-default" "ldb-rand" "ckqt" "moduleCommonLBs" "moduleNDMeshStreamer" "modulecompletion" "hwloc_embedded" "mpi-mainmodule")


#Set up and find boost
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED OFF)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost REQUIRED COMPONENTS system filesystem) 
include_directories(${Boost_INCLUDE_DIRS}) 
link_directories(${Boost_LIBRARY_DIRS}) 

#Set up and find HDF5
set(HDF5_USE_STATIC_LIBS OFF)
find_package(HDF5)
include_directories(${HDF5_INCLUDE_DIRS}) 
add_definitions(${HDF5_DEFINITIONS})

#Set up and find png
find_package(PNG REQUIRED)
include_directories(${PNG_INCLUDE_DIRS}) 
add_definitions(${PNG_DEFINITIONS})

#Set up the linker
set(ENZO_E_LIBRARIES  "dl" "pthread" "z" "gfortran" ${Boost_LIBRARIES}  ${PNG_LIBRARIES} )

#Define preprocessor definitions
add_compile_definitions(CONFIG_PRECISION_DOUBLE)
add_compile_definitions(CONFIG_NODE_SIZE=4)
add_compile_definitions(CONFIG_NODE_SIZE_3=12)
add_compile_definitions(NO_FREETYPE)
add_compile_definitions(SMALL_INTS)
add_compile_definitions(CONFIG_USE_PERFORMANCE CONFIG_USE_MEMORY CONFIG_HAVE_VERSION_CONTROL)

#Define function to process Charm++ interface files
function(preprocess_charm_interfaces out_headers)
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
  set(${out_headers} "${result}" PARENT_SCOPE)
endfunction()

#Copy the input files to the build directory
file(COPY ${CMAKE_SOURCE_DIR}/input DESTINATION ${CMAKE_BINARY_DIR})

#Set some C,CXX, Fortran flags for compiling
#FIXME: Shouldn't need to define these explicitly
set( CMAKE_C_FLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic   ")
set( CMAKE_CXX_FLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic" )
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffast-math -funroll-loops -fPIC -pedantic -ffixed-line-length-132    -fdefault-real-8 -fdefault-double-8")
set( CMAKE_EXE_LINKER_FLAGS "-O3 -g -ffast-math -funroll-loops -fPIC -pedantic -module CommonLBs" )


SET(CMAKE_LINKER ${CHARMC_EXE})
set(CMAKE_CXX_LINK_EXECUTABLE "<CMAKE_LINKER> ${CMAKE_EXE_LINKER_FLAGS} <OBJECTS> <LINK_LIBRARIES> -o <TARGET>")

add_subdirectory(src)

