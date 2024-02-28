# See LICENSE_CELLO file for license and copyright information

# add source files for 

# create Enzo::particle, which includes code involving feedback-related
# particles (e.g. star particles) to the enzo target
#
# We explicitly list files (rather than use glob) since that makes CMake better
# at rebuilds (especially after changing branches)
add_library(Enzo_particle
  particle.hpp
  EnzoMethodPmUpdate.cpp EnzoMethodPmUpdate.hpp
  FofLib.cpp FofLib.hpp

  formation/EnzoBondiHoyleSinkParticle.cpp formation/EnzoBondiHoyleSinkParticle.hpp
  formation/EnzoFluxSinkParticle.cpp formation/EnzoFluxSinkParticle.hpp
  formation/EnzoMethodAccretion.cpp formation/EnzoMethodAccretion.hpp
  formation/EnzoSinkParticle.cpp formation/EnzoSinkParticle.hpp
  formation/EnzoMethodBondiHoyleAccretion.cpp formation/EnzoMethodBondiHoyleAccretion.hpp
  formation/EnzoMethodFluxAccretion.cpp formation/EnzoMethodFluxAccretion.hpp
  formation/EnzoMethodMergeSinks.cpp formation/EnzoMethodMergeSinks.hpp
  formation/EnzoMethodSinkMaker.cpp formation/EnzoMethodSinkMaker.hpp
  formation/EnzoMethodStarMaker.cpp formation/EnzoMethodStarMaker.hpp
  formation/EnzoMethodStarMakerSTARSS.cpp formation/EnzoMethodStarMakerSTARSS.hpp
  formation/EnzoMethodStarMakerStochasticSF.cpp formation/EnzoMethodStarMakerStochasticSF.hpp
  formation/EnzoMethodThresholdAccretion.cpp formation/EnzoMethodThresholdAccretion.hpp

  feedback/EnzoMethodDistributedFeedback.cpp feedback/EnzoMethodDistributedFeedback.hpp
  feedback/EnzoMethodFeedback.cpp feedback/EnzoMethodFeedback.hpp
  feedback/EnzoMethodFeedbackSTARSS.cpp feedback/EnzoMethodFeedbackSTARSS.hpp
)
add_library(Enzo::particle ALIAS Enzo_particle)

# the following can be improved (it links more than is necessary)
target_link_libraries(Enzo_particle PUBLIC enzo ${CELLO_LIBS})
target_include_directories(Enzo_particle PUBLIC ${ROOT_INCLUDE_DIR})
target_link_options(Enzo_particle PRIVATE ${Cello_TARGET_LINK_OPTIONS})
