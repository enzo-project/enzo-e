# See LICENSE_CELLO file for license and copyright information

# Add source files related to problem-initialization to the enzo target
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_initial
  initial.hpp
  EnzoInitialBCenter.cpp EnzoInitialBCenter.hpp
  EnzoInitialBurkertBodenheimer.cpp EnzoInitialBurkertBodenheimer.hpp
  EnzoInitialCloud.cpp EnzoInitialCloud.hpp
  EnzoInitialCollapse.cpp EnzoInitialCollapse.hpp
  EnzoInitialCosmology.cpp EnzoInitialCosmology.hpp
  EnzoInitialImplosion2.cpp EnzoInitialImplosion2.hpp
  EnzoInitialInclinedWave.cpp EnzoInitialInclinedWave.hpp
  EnzoInitialIsolatedGalaxy.cpp EnzoInitialIsolatedGalaxy.hpp
  EnzoInitialSedovArray2.cpp EnzoInitialSedovArray2.hpp
  EnzoInitialSedovArray3.cpp EnzoInitialSedovArray3.hpp
  EnzoInitialSedovRandom.cpp EnzoInitialSedovRandom.hpp
  EnzoInitialShockTube.cpp EnzoInitialShockTube.hpp
  EnzoInitialShuCollapse.cpp EnzoInitialShuCollapse.hpp
  EnzoInitialSoup.cpp EnzoInitialSoup.hpp
  EnzoInitialTurbulence.cpp EnzoInitialTurbulence.hpp
  turboinit.F
  turboinit2d.F
)
add_library(Enzo::initial ALIAS Enzo_initial)

# the following can be improved
target_link_libraries(Enzo_initial PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_initial PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_initial PRIVATE ${Cello_TARGET_LINK_OPTIONS})
