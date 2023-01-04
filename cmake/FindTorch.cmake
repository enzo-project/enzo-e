# See LICENSE_ENZO file for license and copyright information

#.rst:
# FindTorch
# -------
#
# Finds the LibTorch library
#
# This will define the following variables::
#
#   Torch_FOUND    - True if the system has the Grackle library
#
# and the following imported targets::
#
#   Torch::Torch - The Torch library

find_path(Torch_INCLUDE_DIR
  NAMES torch/script.h
  PATHS ${Torch_ROOT}/include
)

# Grackle by default builds both dynamic and static libs.
# We'll use the static one.
find_library(Torch_LIBRARY
  NAMES libtorch.so
  PATHS ${Torch_ROOT}/lib
)

# Ideally we'd also set a version, but I'm still looking for it in the file grackle installed
#set(Grackle_VERSION 0.0.0)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Torch
  FOUND_VAR Torch_FOUND
  REQUIRED_VARS
    Torch_LIBRARY
    Torch_INCLUDE_DIR
#  VERSION_VAR Torch_VERSION
)

if(Torch_FOUND)
  set(Torch_LIBRARIES ${Torch_LIBRARY})
  set(Torch_INCLUDE_DIRS ${Torch_INCLUDE_DIR})
endif()

# Add a target so that we can add torch (and `include` folders) via `target_link_libraries()` directly
if(Torch_FOUND AND NOT TARGET Torch::Torch)
  add_library(Torch::Torch UNKNOWN IMPORTED)
  set_target_properties(Torch::Torch PROPERTIES
    IMPORTED_LOCATION "${Torch_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${Torch_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${Torch_LIBRARY}"
  )
endif()

mark_as_advanced(
  Torch_INCLUDE_DIR
  Torch_LIBRARY
)

# compatibility variables
#set(Grackle_VERSION_STRING ${Grackle_VERSION})
