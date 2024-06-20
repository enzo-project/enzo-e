# Handle external dependencies
# -> locate all existing prebuilt dependencies
# -> in the future, we may also define the recipies for building dependencies
#    that are not pre-built

find_package(Charm REQUIRED)
# Link executables with the charmc wrapper
STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>" "${CHARM_LINKER}"
       CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")


# Need to process Grackle here as both Cello and Enzo-E depend on it
option(USE_GRACKLE "Use Grackle Chemistry" ON)
# don't bother advertising the following option (but users can
# overwrite it if they really want to - e.g. to reduce binary size)
option(GRACKLE_USE_STATIC_LIBS  "sets Grackle's lib-type if USE_GRACKLE=ON" ON)
if (USE_GRACKLE)
  find_package(Grackle)
  if (Grackle_FOUND)
    # This really only needs to be defined for a relatively small subset of
    # files in the Enzo-layer
    add_compile_definitions(CONFIG_USE_GRACKLE)

  else()
    message(FATAL_ERROR
      "Configured to use Grackle but Grackle library not found.\n"
      "Either disable grackle (e.g., `-DUSE_GRACKLE=OFF`) or provide path "
      "(e.g., `-DGrackle_ROOT=/PATH/TO/GRACKLE/INSTALL/DIRECTORY`).")
  endif()
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
