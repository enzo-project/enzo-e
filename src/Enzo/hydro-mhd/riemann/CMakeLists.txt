# See LICENSE_CELLO file for license and copyright information

# define the Enzo:riemann subcomponent
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_riemann
  EnzoRiemann.hpp EnzoRiemann.cpp
  EnzoRiemannLUT.hpp
  EnzoRiemannUtils.hpp
  EnzoRiemannHLL.hpp
  EnzoRiemannHLLC.hpp
  EnzoRiemannHLLD.hpp
  EnzoRiemannImpl.hpp
)
add_library(Enzo::riemann ALIAS Enzo_riemann)

# the following links a LOT of unnecessary stuff!
# - in the future, it might be nice to only depend on subcomponents of the
#   Cello or Enzo layer (to reduce how often this library gets recompiled)
target_link_libraries(Enzo_riemann PUBLIC enzo ${Cello_LIBS})
target_include_directories (Enzo_riemann PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_riemann PRIVATE ${Cello_TARGET_LINK_OPTIONS})
