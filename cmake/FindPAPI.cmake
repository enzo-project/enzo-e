# See LICENSE_ENZO file for license and copyright information

#.rst:
# FindPAPI
# -------
#
# Finds the PAPI performance API
#
# This will define the following variables::
#
#   PAPI_FOUND    - True if the system has the PAPI library
#
# and the following imported targets::
#
#   PAPI::papi - The PAPI library

find_path(PAPI_INCLUDE_DIR
  NAMES papi.h
)

find_library(PAPI_LIBRARY
  NAMES papi
)

# Ideally we'd also set a version, but I'm still looking for it in the file papi installed
#set(PAPI_VERSION 0.0.0)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PAPI
  FOUND_VAR PAPI_FOUND
  REQUIRED_VARS
    PAPI_LIBRARY
    PAPI_INCLUDE_DIR
#  VERSION_VAR PAPI_VERSION
)

if(PAPI_FOUND)
  set(PAPI_LIBRARIES ${PAPI_LIBRARY})
  set(PAPI_INCLUDE_DIRS ${PAPI_INCLUDE_DIR})
endif()

# Add a target so that we can add papi (and `include` folders) via `target_link_libraries()` directly
if(PAPI_FOUND AND NOT TARGET PAPI::papi)
  add_library(PAPI::papi UNKNOWN IMPORTED)
  set_target_properties(PAPI::papi PROPERTIES
    IMPORTED_LOCATION "${PAPI_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${PAPI_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${PAPI_LIBRARY}"
  )
endif()

mark_as_advanced(
  PAPI_INCLUDE_DIR
  PAPI_LIBRARY
)

# compatibility variables
#set(PAPI_VERSION_STRING ${PAPI_VERSION})
