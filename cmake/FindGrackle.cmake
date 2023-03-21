# See LICENSE_ENZO file for license and copyright information

#.rst:
# FindGrackle
# -------
#
# Finds the Grackle chemistry library
#
# Results
# =======
#
# This will define the following variable:
#
#   Grackle_FOUND    - True if the system has the Grackle library
#
# and the following imported targets:
#
#   Grackle::Grackle - The Grackle library target
#
# If possible, this module will also determine the Grackle version. The
# functionality to query this information was not provided in any stable
# releases of Grackle before version 3.2, so it probably won't be reported for
# earlier versions.
#
#   Grackle_VERSION  - Encodes the library version number.
#
# Note that the version numbers for Grackle's development-versions (which users
# are encouraged to use) don't strictly follow CMake's expectations about
# version numbers. Thus, developers should probably refrain from using
# ``find_package``'s ability to check version compatability.
#
# Also note that searching for the Grackle version requires a call to
# ``try_run``, which can produce challenges during cross-compilation. For that
# reason, ``GRACKLE_SKIP_VERSION_SEARCH`` can be set to ``ON`` to avoid
# searching for the version number.
#
# Advanced Usage
# ==============
# Set GRACKLE_USE_STATIC_LIBS to ``ON`` to indicate preference for static
# libraries. The default is ``OFF``.

function(_Grackle_identify_version_and_store)
  # This function searches for the version number. Upon success, it will modify
  # the Grackle_VERSION variable in the main scope.
  #
  # This is primarily written as a function to simplify internal variable names
  set(test_file
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/grackle/get_vers.c)
  file(WRITE ${test_file}
    "#ifdef __cplusplus\n"
    "#include <cstdio>\n"
    "extern \"C\" {\n"
    "#include <grackle.h>\n"
    "}\n"
    "#else\n"
    "#include <stdio.h>\n"
    "#include <grackle.h>\n"
    "#endif\n"
    "int main(int argc, char* argv[]){\n"
    "  grackle_version gversion = get_grackle_version();\n"
    "#ifdef __cplusplus\n"
    "  std::puts(gversion.version);\n"
    "#else\n"
    "  puts(gversion.version);\n"
    "#endif\n"
    "  return 0;\n"
    "}")

  try_run(exit_code compile_success ${CMAKE_BINARY_DIR} ${test_file}
    LINK_LIBRARIES Grackle::Grackle
    RUN_OUTPUT_VARIABLE test_output
  )

  if(compile_success) # successful compilation

    if (exit_code EQUAL 0) # run succeeded
      string(STRIP "${test_output}" striped_version_str)
      set(Grackle_VERSION "${striped_version_str}" PARENT_SCOPE)
    else ()
      if ("${exit_code}" STREQUAL "FAILED_TO_RUN")
        set(err_explanation "failed to run")
      else ()
        set(err_explanation "had the exit code ${exit_code}")
      endif()
      message(FATAL_ERROR
        "Something unexpected happened while trying to run the test program "
        "for querying Grackle's version. The test program compiled, but "
        "the test program ${err_explanation}. You can avoid searching for "
        "Grackle's version number by setting GRACKLE_SKIP_VERSION_SEARCH=ON")
    endif()
  endif()
endfunction()

find_path(Grackle_INCLUDE_DIR
  NAMES grackle.h
)

if(DEFINED CMAKE_FIND_LIBRARY_SUFFIXES)
  set(_grackle_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES
      "${CMAKE_FIND_LIBRARY_SUFFIXES}")
else()
  set(_grackle_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES)
endif()

# adjust CMAKE_FIND_LIBRARY_SUFFIXES if there’s a preference for static
# libraries
if(GRACKLE_USE_STATIC_LIBS)
  if(WIN32) # this branch isn't relevant, but it's part of the standard idiom
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
  endif()
endif()

find_library(Grackle_LIBRARY
  NAMES grackle
)

# restore initial value of CMAKE_FIND_LIBRARY_SUFFIXES
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_grackle_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})


if (Grackle_LIBRARY AND Grackle_INCLUDE_DIR)

  # construct imported Grackle target BEFORE find_package_handle_standard_args.
  # - There's precedent for doing this in FindPNG.cmake and FindMPI.cmake
  #   modules that are shipped with CMake
  # - This is necessary if we want to query Grackle's version
  if(NOT TARGET Grackle::Grackle)
    add_library(Grackle::Grackle UNKNOWN IMPORTED)
    set_target_properties(Grackle::Grackle PROPERTIES
      IMPORTED_LOCATION "${Grackle_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${Grackle_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES "${Grackle_LIBRARY}"
    )
  endif()

  # try to determine the version
  message("${CMAKE_FILES_DIRECTORY}")
  if (NOT GRACKLE_SKIP_VERSION_SEARCH)
    _Grackle_identify_version_and_store()
  endif()

endif()

# if we can't find a version number (or we choose not to look), the
# Grackle_VERSION variable will be undefined. It seems to be ok to pass
# Grackle_VERSION (when it's undefined) to find_package_handle_standard_args;
# the desired behavior is achieved (FindMPI.cmake also does the same thing)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Grackle
  FOUND_VAR Grackle_FOUND
  REQUIRED_VARS
    Grackle_LIBRARY
    Grackle_INCLUDE_DIR
  VERSION_VAR Grackle_VERSION
)

if(Grackle_FOUND)
  set(Grackle_LIBRARIES ${Grackle_LIBRARY})
  set(Grackle_INCLUDE_DIRS ${Grackle_INCLUDE_DIR})
endif()

mark_as_advanced(
  Grackle_INCLUDE_DIR
  Grackle_LIBRARY
)
