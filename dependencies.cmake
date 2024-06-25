# Handle external dependencies
# -> locate all existing prebuilt dependencies
# -> in the future, we may also define the recipies for building dependencies
#    that are not pre-built

find_package(Charm REQUIRED)
# Link executables with the charmc wrapper
STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>" "${CHARM_LINKER}"
       CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")

# ---------------------------------------------------------------------------
# Next, handle all dependencies that may be compiled as part of this build
# ---------------------------------------------------------------------------
# If we want to build these dependencies, there are 2 ways to do this:
# 1. assume that the source code is embedded within the GitHub repository in
#    a known location (this is usually accomplished via git submodules)
# 2. automatically download the source-code of the dependency with cmake's
#    FetchContent functionality
#
# For now, we adopt the FetchContent functionality since the cmake devs view
# that as the "canonical" way to do things
# -> once we can require that everyone uses CMake 3.24.2 or newer, we can
#    leverage some useful features built into FetchContent that will make
#    swapping between find_package and FetchContent much more seamless
#
# NOTE: when using FetchContent with git repositories, it is considered
# idiomatic to use commit-hashes rather than tag-names or commit names.
# -> the reason is that hashes can't change (while the other 2 choices can)
# -> if we don't use a hash, then cmake will contact the remote repository
#    between every build to check whether the tag or branch name

include(FetchContent)

# Need to process Grackle here as both Cello and Enzo-E depend on it
option(USE_GRACKLE "Use Grackle Chemistry" ON)
option(USE_EXTERNAL_GRACKLE "Use external grackle. Meaningless unless USE_GRACKLE is ON" OFF)
if(USE_GRACKLE)
  # first, declare where to find FetchContent would find grackle (even if we
  # don't ultimately use FetchContent

  # TODO: shift this to mainline Grackle Repository after
  #         https://github.com/grackle-project/grackle/pull/182
  #       get's merged into Grackle
  FetchContent_Declare(Grackle
    GIT_REPOSITORY https://github.com/mabruzzo/grackle
    GIT_TAG 689be185ac55dba098309e2da9d6acdda37d1923
    # ^ current hash is right after cmake build-system got introduced
  )

  # some time after the following PR is merged
  #    https://github.com/grackle-project/grackle/pull/204
  # And we feel comfortable requiring that people use a CMake-build of grackle
  # that ships a GrackleConfig.cmake file, a better default behavior might be
  # to always search for Grackle first and then perform an in-source build if
  # it can't be found

  if (USE_EXTERNAL_GRACKLE)
    find_package(Grackle)

    if (NOT Grackle_FOUND)
      message(FATAL_ERROR
        "Configured to use Grackle but Grackle library not found.\n"
        "Either disable grackle (e.g. `-DUSE_GRACKLE=OFF`), start using "
        "automatic dependency management (e.g. `-DUSE_EXTERNAL_GRACKLE=OFF), "
        "or provide path "
        "(e.g. `-DGrackle_ROOT=/PATH/TO/GRACKLE/INSTALL/DIRECTORY`).")
    endif()
  else()
    set(GRACKLE_USE_DOUBLE USE_DOUBLE_PREC)
    set(GRACKLE_USE_OPENMP OFF)

    message("downloading Grackle...")
    FetchContent_MakeAvailable(Grackle)
  endif()

  # This really only needs to be defined for a relatively small subset of
  # files in the Enzo-layer
  add_compile_definitions(CONFIG_USE_GRACKLE)
endif()


find_package(PNG REQUIRED)
add_compile_definitions(NO_FREETYPE)

find_package(HDF5 REQUIRED COMPONENTS C)
# HDF5 Interface library
add_library(HDF5_C INTERFACE)
target_link_libraries(HDF5_C INTERFACE ${HDF5_C_LIBRARIES})
target_compile_definitions(HDF5_C INTERFACE ${HDF5_C_DEFINITIONS})
target_include_directories(HDF5_C INTERFACE ${HDF5_C_INCLUDE_DIRS})

option (use_jemalloc "Use the jemalloc library for memory allocation" OFF)
if (use_jemalloc)
  find_package(jemalloc)
  if (jemalloc_FOUND)
    add_compile_definitions(CONFIG_USE_JEMALLOC)
  else()
    message(FATAL_ERROR
      "Requested to use the jemalloc library for memory allocation but jemalloc was not found. "
      "Try setting specific path via `-Djemalloc_ROOT=/PATH/TO/jemalloc/INSTALL` "
      " or disable jemalloc via `-Duse_jemalloc=OFF` (default)."
      )
  endif()
endif()

option(use_papi "Use the PAPI performance API" OFF)
if (use_papi)
  find_package(PAPI)
  if (PAPI_FOUND)
    add_compile_definitions(CONFIG_USE_PAPI PAPI3)
  else()
    message(FATAL_ERROR
      "Requested to use PAPI performance API but PAPI was not found. "
      "Try setting specific path via `-DPAPI_ROOT=/PATH/TO/PAPI/INSTALL` "
      " or disable PAPI via `-Duse_papi=OFF` (default)."
      )
  endif()
endif()
