# See LICENSE_ENZO file for license and copyright information

#.rst:
# Findjemalloc
# -------
#
# Finds the jemalloc for memory allocation
#
# This will define the following variables::
#
#   jemalloc_FOUND    - True if the system has the jemalloc library
#
# and the following imported targets::
#
#   jemalloc::jemalloc - The jemalloc library

find_path(jemalloc_INCLUDE_DIR
  NAMES jemalloc.h jemalloc/jemalloc.h
)

find_library(jemalloc_LIBRARY
  NAMES jemalloc
)

# Ideally we'd also set a version, but I'm still looking for it in the file jemalloc installed
#set(jemalloc_VERSION 0.0.0)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(jemalloc
  FOUND_VAR jemalloc_FOUND
  REQUIRED_VARS
    jemalloc_LIBRARY
    jemalloc_INCLUDE_DIR
#  VERSION_VAR jemalloc_VERSION
)

if(jemalloc_FOUND)
  set(jemalloc_LIBRARIES ${jemalloc_LIBRARY})
  set(jemalloc_INCLUDE_DIRS ${jemalloc_INCLUDE_DIR})
endif()

# Add a target so that we can add jemalloc (and `include` folders) via `target_link_libraries()` directly
if(jemalloc_FOUND AND NOT TARGET jemalloc::jemalloc)
  add_library(jemalloc::jemalloc UNKNOWN IMPORTED)
  set_target_properties(jemalloc::jemalloc PROPERTIES
    IMPORTED_LOCATION "${jemalloc_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${jemalloc_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${jemalloc_LIBRARY}"
  )
endif()

mark_as_advanced(
  jemalloc_INCLUDE_DIR
  jemalloc_LIBRARY
)

# compatibility variables
#set(jemalloc_VERSION_STRING ${jemalloc_VERSION})
