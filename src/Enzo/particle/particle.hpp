// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/particle/particle.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for particle subcomponent within the \ref Enzo layer

#ifndef ENZO_PARTICLE_PARTICLE_HPP
#define ENZO_PARTICLE_PARTICLE_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>
#include <vector>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Method

#include "Enzo/enzo.hpp" // enzo_float, EnzoBlock

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

#include "particle/EnzoMethodPmUpdate.hpp"

// [order dependencies:]
#include "particle/formation/EnzoSinkParticle.hpp"
#include "particle/formation/EnzoBondiHoyleSinkParticle.hpp"
#include "particle/formation/EnzoFluxSinkParticle.hpp"

#include "particle/formation/EnzoMethodAccretion.hpp"
#include "particle/formation/EnzoMethodBondiHoyleAccretion.hpp"
#include "particle/formation/EnzoMethodFluxAccretion.hpp"
#include "particle/formation/EnzoMethodMergeSinks.hpp"
#include "particle/formation/EnzoMethodSinkMaker.hpp"
#include "particle/formation/EnzoMethodStarMaker.hpp"
#include "particle/formation/EnzoMethodStarMakerSTARSS.hpp"
#include "particle/formation/EnzoMethodStarMakerStochasticSF.hpp"
#include "particle/formation/EnzoMethodThresholdAccretion.hpp"

#include "particle/feedback/EnzoMethodDistributedFeedback.hpp"
#include "particle/feedback/EnzoMethodFeedback.hpp"
#include "particle/feedback/EnzoMethodFeedbackSTARSS.hpp"

#endif /* ENZO_PARTICLE_PARTICLE_HPP */
