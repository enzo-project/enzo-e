# See LICENSE_ENZO file for license and copyright information

#.rst:
# CreateStdFilesystemTarget
# -------------------------
#
# This module defines the create_StdFilesystem_target function
#
# That function:
# - confirms that the chosen compiler supports the file system library included
#   in the standard library starting in C++17
# - defines the StdFilesystem::StdFilesystem library to communicate extra usage
#   requirements
#
# In 99% of cases, there shouldn't be any usage requirements. The main
# exception relates to the implementations of the standard libraries that were
# available around the time that the c++17 standard was ratified. For example:
# - libstdc++ (GNU's implementation) required the `-lstdc++fs` flag
# - libc++ (LLVM's implementation) required the `-lc++fs` flag
# This extra argument mirrors the way in which implementations of libc on unix
# systems may require the `-lm` flag when using math functions.
#
# At the time of writing this, much of the extra functionality for enabling
# compatability with these older versions is hypothetical (it hasn't been
# tested).
#
# FUTURE WORK
# It may be worth asking "do we want to support these older versions?" Or do we
# simply want to test if things "just work" out of the box? If things break in
# the latter case, we could just say what's probably going wrong and tell the
# user that they probably need to upgrade. As time progresses, the latter
# option will become more appealing (and delete ~50% of this file's logic)

function(_detect_stdlib_impl output_variable)
  # writes the name of the standard library implementation to the variable name
  # specifed by output_variable
  # -> known options include
  #    "libstdc++" (The GNU implementation)
  #    "libc++" (The LLVM implementation)
  #    "unknown"
  # -> Note: the standard library doesn't have to match the compiler (e.g.
  #    clang on linux often uses libstdc++)

  set(test_file
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/std-fs/stdimpl.cpp)

  # this detection stategy, won't work for libstdc++ older than 6.1
  # (but that's okay, such versions are older than C++17)
  file(WRITE ${test_file}
" 
// include headers defining implementation-specific macros (since we include
// <iostream>, this may not be strictly necessary)
#if __has_include(<ciso646>)
  // pre c++20, this header is commonly included for std library implementation 
  // macros: https://en.cppreference.com/w/cpp/header/version
  #include <ciso646>
#else
  #include <version>
#endif

#include <iostream>
int main(int argc, char* argv[]){
#ifdef __GLIBCXX__
  std::cout << \"libstdc++\\n\";
#elif defined(_LIBCPP_VERSION)
  std::cout << \"libc++\\n\";
#else
  std::cout << \"unknown\\n\";
#endif
  return 0;
}")

  try_run(exit_code compile_success ${CMAKE_BINARY_DIR} ${test_file}
    RUN_OUTPUT_VARIABLE test_output
  )

  if(compile_success AND (exit_code EQUAL 0)) # successful test
    string(STRIP "${test_output}" stripped_implementation_name)
    set(${output_variable} "${stripped_implementation_name}" PARENT_SCOPE)
  else()
    message(FATAL_ERROR "Something unexpected happened")
  endif()
endfunction()


function(_StdFilesystem_try_compile link_library output_variable)
  # compile & run a test program with std::filesystem
  # -> stores a value of 1 or 0 (denoting success) in the output_variable

  set(test_file
    ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/std-fs/test.cpp)

  # the whole point is for this to be extremely simple
  file(WRITE ${test_file}
"
#include <filesystem>
#include <iostream>

int main(int argc, char* argv[]){
  std::cout << std::filesystem::temp_directory_path();
  return 0;
}")

  if (link_library STREQUAL "")
    try_run(exit_code compile_success ${CMAKE_BINARY_DIR} ${test_file}
      RUN_OUTPUT_VARIABLE test_output
    )
  else()
    try_run(exit_code compile_success ${CMAKE_BINARY_DIR} ${test_file}
      LINK_LIBRARIES ${link_library}
      RUN_OUTPUT_VARIABLE test_output
    )
  endif()

  if(compile_success AND (exit_code EQUAL 0)) # successful test
    set(${output_variable} "1" PARENT_SCOPE)
  else()
    set(${output_variable} "0" PARENT_SCOPE)
  endif()
endfunction()

function(_create_StdFilesystem_target_helper link_library)
  add_library(StdFilesystem::StdFilesystem INTERFACE IMPORTED)
  if(NOT ${link_library} STREQUAL "")
    target_link_libraries(StdFilesystem::StdFilesystem
      INTERFACE "${link_library}")
  endif()
endfunction()

function(create_StdFilesystem_target)
  # this does the heavy lifting
  if(TARGET StdFilesystem::StdFilesystem)
    return()
  elseif(DEFINED CACHE{__StdFilesystem_LINKLIBRARY})
    _create_StdFilesystem_target_helper("${__StdFilesystem_LINKLIBRARY}")
    return()
  endif()

  set(BASE_MSG "Checking std::filesystem support")
  message(STATUS "${BASE_MSG}")

  # first, we try without any link-library
  set(link_library "")
  set(compile_success "0")
  _StdFilesystem_try_compile("${link_library}" compile_success)

  if (compile_success)
    set(SUCCESS_MSG "works out of the box")
  else()
    # first provide a nice error message for specific compilers
    if(CMAKE_CXX_COMPILER_VERSION STREQUAL "AppleClang")
      if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "11.0.0")
        message(FATAL_ERROR
          "Version ${CMAKE_CXX_COMPILER_VERSION} of Apple's Clang compiler "
          "does NOT support std::filesystem. Please upgrade to version "
          "11.0.0 or newer (you probably have to update Xcode). Alternatively "
          "install a different c++ compiler on your system")
      else()
        message(FATAL_ERROR "Unclear why Apple's Clang compiler isn't working")
      endif()
    endif()

    # gracefully handle cases where the user has an older version of a
    # standard-library implementation that requires extra compiler args
    # - NOTE: the precise implementation does not necessarily correspond to the
    #   compiler (e.g. GNU's implementation is commonly used with clang++)
    _detect_stdlib_impl(stdlib_impl)

    if (stdlib_impl STREQUAL "libstdc++")
      # this is the GNU implementation of the standard library
      set(link_library "stdc++fs")
      set(min_req_vers "8.0.0")
      set(link_lib_unneeded_vers "9.1")
      _StdFilesystem_try_compile("${link_library}" compile_success)
    elseif (stdlib_impl STREQUAL "libc++")
      # this is the LLVM implementation of the standard library
      set(link_library "c++fs")
      set(min_req_vers "7")
      set(link_lib_unneeded_vers "9")
      _StdFilesystem_try_compile("${link_library}" compile_success)
    else()
      message(FATAL_ERROR
        "std::filesystem doesn't work and the implementation of the standard"
        "library can't be determined (it's not a recent version of libstdc++ "
        "or libc++). It's likely that: the implementation doesn't support "
        "std::filesystem OR we may need to pass an extra linker flag to use "
        "std::filesystem (as is required by earlier versions of libstdc++ or "
        "libc++)")
    endif()

    if(NOT compile_success)
      message(FATAL_ERROR
        "std::filesystem doesn't work. You appear to be using the "
        "${stdlib_impl} implementation of the standard library. Note that a "
        "version of at least ${min_req_vers} is required.")
    endif()
    set(SUCCESS_MSG
      "works with the linker flag -l${link_library} (this is requirement of "
      "${stdlib_impl} before version ${link_lib_unneeded_vers}")
  endif()

  message(STATUS "${BASE_MSG} -- ${SUCCESS_MSG}")

  set(__StdFilesystem_LINKLIBRARY "${link_library}" CACHE INTERNAL
    "empty-string or library we must link against to use std::filesystem")
  _create_StdFilesystem_target_helper("${__StdFilesystem_LINKLIBRARY}")

endfunction()
