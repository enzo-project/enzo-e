// See LICENSE_CELLO file for license and copyright information

/// @file     Enzo/fluid-props/fluid-props.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-05-30
/// @brief    Include file for fluid-props subcomponent of the \ref Enzo layer

#ifndef ENZO_FLUIDPROPS_FLUIDPROPS_HPP
#define ENZO_FLUIDPROPS_FLUIDPROPS_HPP

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <string>

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "Cello/cello.hpp"

#include "Cello/compute.hpp" // Compute
#include "Cello/mesh.hpp"    // Block
#include "Cello/problem.hpp" // Physics
#include "Cello/view.hpp" // CelloView

#include "Enzo/enzo.hpp" // enzo_float, EnzoBlock, EnzoEFltArrayMap

//----------------------------------------------------------------------
// Component headers
//----------------------------------------------------------------------

// [order dependencies:]
#include "fluid-props/EnzoEOSIdeal.hpp"
#include "fluid-props/EnzoEOSIsothermal.hpp"
#include "fluid-props/EnzoEOSVariant.hpp"

#include "fluid-props/EnzoDualEnergyConfig.hpp"
#include "fluid-props/EnzoFluidFloorConfig.hpp"
#include "fluid-props/EnzoPhysicsFluidProps.hpp"

#include "fluid-props/EnzoComputePressure.hpp"
#include "fluid-props/EnzoComputeTemperature.hpp"

#endif /* ENZO_FLUIDPROPS_FLUIDPROPS_HPP */
