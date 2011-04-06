// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_HPP
#define ENZO_HPP

/// @file     enzo.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Include file for the \ref Enzo component

//----------------------------------------------------------------------
// System includes
//----------------------------------------------------------------------

#include <vector>
#include <string>
#include <limits>

//----------------------------------------------------------------------
// Cello include file
//----------------------------------------------------------------------

#include "cello.hpp"

//----------------------------------------------------------------------
// Component dependencies
//----------------------------------------------------------------------

#include "performance.hpp"
#include "method.hpp"
#include "simulation.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "enzo_defines.hpp"
#include "enzo_typedefs.hpp"
#include "enzo_fortran.hpp"

#include "enzo_EnzoNamespace.hpp"

#include "enzo_EnzoFactory.hpp"

#include "enzo_EnzoSimulation.hpp"
#   include "enzo_EnzoSimulationMpi.hpp"
#   include "enzo_EnzoSimulationCharm.hpp"

#include "enzo_EnzoBlock.hpp"

#include "enzo_EnzoTimestep.hpp"
#include "enzo_EnzoBoundary.hpp"
#include "enzo_EnzoInitialImplosion2.hpp"
#include "enzo_EnzoMethodPpm.hpp"
#include "enzo_EnzoMethodPpml.hpp"
#include "enzo_EnzoOutputImage.hpp"

#endif /* ENZO_HPP */

