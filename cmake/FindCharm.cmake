################################################################################
#
# \file      FindCharm.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved.
#
# This program was produced under U.S. Government contract 89233218CNA000001 for
# Los Alamos National Laboratory (LANL), which is operated by Triad National
# Security, LLC for the U.S. Department of Energy/National Nuclear Security
# Administration. All rights in the program are reserved by Triad National
# Security, LLC, and the U.S. Department of Energy/National Nuclear Security
# Administration. The Government is granted for itself and others acting on its
# behalf a nonexclusive, paid-up, irrevocable worldwide license in this material
# to reproduce, prepare derivative works, distribute copies to the public,
# perform publicly and display publicly, and to permit others to do so.
# 
# Additionally, redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation and/or
# other materials provided with the distribution.
# 
# 3. Neither the name of Triad National Security, LLC, Los Alamos National
# Laboratory, LANL, the U.S. Government, nor the names of its contributors may be
# used to endorse or promote products derived from this software without specific
# prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY TRIAD NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL SECURITY, LLC OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
# GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# \brief     Find Charm++ from https://github.com/quinoacomputing/cmake-modules
#
# 2022/07/18 mabruzzo: refactoring logic for detecting whether Charm++ was compiled with SMP mode
#
# 2023/02/05 mabruzzo: add logic to export a target for including charm++ headers
################################################################################

# Charm++: http://charmplusplus.org
#
#  CHARM_FOUND          - True if the charmc compiler wrapper was found
#  CHARM_INCLUDE_DIRS   - Charm++ include files paths
#
#  CHARM::CHARM_HEADERS - target mostly just used for book-keeping purposes.
#                         This manages include directories. It doesn't actually
#                         specify library locations. Instead, that info is
#                         automatically provided when using CHARM_COMPILER for
#                         linking (that is shown down below)
#
#  Set CHARM_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(CHARM_ROOT "/path/to/custom/charm") # prefer over system
#  find_package(Charm)
#  if(CHARM_FOUND)
#    # Link executables with the charmc wrapper
#    STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>" "${CHARM_COMPILER}"
#           CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")
#  endif()

function(_GET_CHARMINC _OUT_INC _charmc)
  file(STRINGS ${_charmc} _contents REGEX "^CHARMINC=")
  if(_contents)
    string(REGEX REPLACE "^CHARMINC=\"(.*)\"" "\\1" ${_OUT_INC} "${_contents}")
    set(${_OUT_INC} ${${_OUT_INC}} PARENT_SCOPE)
  else()
    message(FATAL_ERROR "file ${_charmc} does not exist")
  endif()
endfunction()

# find out if Charm++ was built with randomzied message queues
function(GET_CHARM_QUEUE_TYPE CHARM_RNDQ conv_conf_hdr)
  file( STRINGS ${conv_conf_hdr} _contents
        REGEX ".*#define CMK_RANDOMIZED_MSGQ[ \t]+" )
  if(_contents)
    string(REGEX REPLACE ".*#define CMK_RANDOMIZED_MSGQ[ \t]+([01]+).*" "\\1" RNDQ "${_contents}")
    if (RNDQ EQUAL 1)
      set(CHARM_RNDQ true PARENT_SCOPE)
      message(STATUS "Charm++ built with randomized message queues")
    endif()
  else()
    message(FATAL_ERROR "Include file ${conv_conf_hdr} does not exist")
 endif()
endfunction()

# If already in cache, be silent
if (CHARM_INCLUDE_DIRS AND CHARM_COMPILER AND CHARM_RUN)
  set (CHARM_FIND_QUIETLY TRUE)
endif()

INCLUDE(FindCygwin)

FIND_PROGRAM(CHARM_COMPILER
  NAMES charmc
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

FIND_PROGRAM(CHARM_RUN
  NAMES charmrun
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

FIND_PROGRAM(AMPI_C_COMPILER
  NAMES ampicc
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

FIND_PROGRAM(AMPI_CXX_COMPILER
  NAMES ampicxx
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

FIND_PROGRAM(AMPI_RUN
  NAMES ampirun
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

if(CHARM_COMPILER)
  _GET_CHARMINC(HINTS_CHARMINC ${CHARM_COMPILER})
endif()

FIND_PATH(CHARM_INCLUDE_DIR NAMES charm.h
                            HINTS ${HINTS_CHARMINC}
                                  ${CHARM_ROOT}/include
                                  $ENV{CHARM_ROOT}/include
                                  ${CMAKE_INSTALL_PREFIX}/charm/include
                            PATH_SUFFIXES charm)

if(CHARM_INCLUDE_DIR)
  FIND_PATH(CHARM_CONV_HDR NAMES conv-autoconfig.h
                           HINTS ${HINTS_CHARMINC}
                                 ${CHARM_INCLUDE_DIR}
                           PATH_SUFFIXES ../tmp)
  GET_CHARM_QUEUE_TYPE(CHARM_QUEUE_TYPE ${CHARM_CONV_HDR}/conv-autoconfig.h)
endif()

if(CHARM_INCLUDE_DIR)
  set(CHARM_INCLUDE_DIRS ${CHARM_INCLUDE_DIR})
else()
  set(CHARM_INCLUDE_DIRS "")
endif()

# Handle the QUIETLY and REQUIRED arguments and set CHARM_FOUND to TRUE if all
# listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Charm DEFAULT_MSG CHARM_COMPILER
                                  CHARM_INCLUDE_DIRS CHARM_RUN)

if(AMPI_C_COMPILER AND AMPI_CXX_COMPILER)
  set(AMPI_FOUND true)
  message(STATUS "Charm++ built with AMPI")
endif()


include(CheckIncludeFiles)
include(CheckCXXSourceCompiles)

if(CHARM_COMPILER AND NOT DEFINED CHARM_SMP)
  # Check whether Charm++ was built in SMP mode. Restructuring of Charm++ (that
  # seems to have occured around the release of charm++ 7.0) makes this test
  # somewhat non-trivial
  #
  # - in earlier versions, the CMK_SMP macro is only defined in the
  #   "conv-mach-opt.h" header if charm++ was built in SMP mode (we could use
  #   CHECK_SYMBOL_EXISTS to query if it was defined) Note: when the macro
  #   is defined, it has a value of 1.
  # - in modern versions, the CMK_SMP macro is ALWAYS defined in the
  #   "conv-autoconfig.h" header. When built in SMP mode, it has a value of 1.
  #   Otherwise, it has a value of 0.

  # construct the include statements for the test program
  set(CHARM_SMP_INCLUDE_STATEMENTS "")
  foreach(header in ITEMS "conv-mach-opt.h" "conv-autoconfig.h")
    CHECK_INCLUDE_FILES("${CHARM_INCLUDE_DIR}/${header}"
                        HAVE_CHARM_SMP_HEADER)
    if (HAVE_CHARM_SMP_HEADER)
      string(APPEND CHARM_SMP_INCLUDE_STATEMENTS
             "#include \"${CHARM_INCLUDE_DIR}/${header}\"\n")
    endif()
    unset(HAVE_CHARM_SMP_HEADER CACHE)
  endforeach()

  # run the test program to determine whether SMP mode is used
  CHECK_CXX_SOURCE_COMPILES("
    ${CHARM_SMP_INCLUDE_STATEMENTS}

    #ifdef CMK_SMP
    #if CMK_SMP > 0
    int main(void) { return 0; }
    #endif
    #endif
    " CHARM_SMP)

  if (CHARM_SMP)
    message(STATUS "Charm++ built in SMP mode")
  else()
    message(STATUS "Charm++ built in non-SMP mode")
  endif()

endif()

# define a target primarily for book-keeping purposes
# - this target provides appropriate include-directories used in linking
# - it doesn't currently provide any information about the library location
#   (this file currently assumes that ${CHARM_COMPILER} is used for linking,
#    which provides that information automatically)
if(CHARM_COMPILER)
  if (NOT TARGET CHARM::CHARM_HEADERS)

    add_library(CHARM::CHARM_HEADERS INTERFACE IMPORTED)
    target_include_directories(CHARM::CHARM_HEADERS
      INTERFACE ${CHARM_INCLUDE_DIRS}
      )

  endif()

endif()

MARK_AS_ADVANCED(CHARM_COMPILER CHARM_INCLUDE_DIRS CHARM_RUN AMPI_FOUND
                 AMPI_C_COMPILER AMPI_CXX_COMPILER AMPI_RUN CHARM_SMP
                 CHARM_QUEUE_TYPE)
