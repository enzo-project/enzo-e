# See LICENSE_ENZO file for license and copyright information

#.rst:
# FindGrackle
# -------
#
# Finds the Grackle chemistry library
#
# This will define the following variables::
#
#   Grackle_FOUND    - True if the system has the Grackle library
#
# and the following imported targets::
#
#   Grackle::Grackle - The Grackle library

find_path(Grackle_INCLUDE_DIR
  NAMES grackle.h
)

# Grackle by default builds both dynamic and static libs.
# We'll use the static one.
find_library(Grackle_LIBRARY
  NAMES grackle
)

# Ideally we'd also set a version, but I'm still looking for it in the file grackle installed
#set(Grackle_VERSION 0.0.0)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Grackle
  FOUND_VAR Grackle_FOUND
  REQUIRED_VARS
    Grackle_LIBRARY
    Grackle_INCLUDE_DIR
#  VERSION_VAR Grackle_VERSION
)

if(Grackle_FOUND)
  set(Grackle_LIBRARIES ${Grackle_LIBRARY})
  set(Grackle_INCLUDE_DIRS ${Grackle_INCLUDE_DIR})
endif()

# Add a target so that we can add grackle (and `include` folders) via `target_link_libraries()` directly
if(Grackle_FOUND AND NOT TARGET Grackle::Grackle)
  add_library(Grackle::Grackle UNKNOWN IMPORTED)
  set_target_properties(Grackle::Grackle PROPERTIES
    IMPORTED_LOCATION "${Grackle_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${Grackle_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${Grackle_LIBRARY}"
  )
endif()

mark_as_advanced(
  Grackle_INCLUDE_DIR
  Grackle_LIBRARY
)

# compatibility variables
#set(Grackle_VERSION_STRING ${Grackle_VERSION})
