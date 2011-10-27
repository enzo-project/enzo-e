// See LICENSE_CELLO file for license and copyright information

/// @file     enzo.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Include file for the \ref Enzo component

#ifndef ENZO_HPP
#define ENZO_HPP

//----------------------------------------------------------------------

#include "fortran.h" /* included so scons knowns to install fortran.h */
#include "enzo_defines.hpp"
#include "enzo_typedefs.hpp"
#include "enzo_fortran.hpp"

//----------------------------------------------------------------------

#define OMEGA_TOLERANCE 1.0e-5
 
#ifdef CONFIG_PRECISION_SINGLE
#   define ETA_TOLERANCE 1.0e-5
#endif
#ifdef CONFIG_PRECISION_DOUBLE
#   define ETA_TOLERANCE 1.0e-10
#endif
#ifdef CONFIG_PRECISION_QUAD
#   define ETA_TOLERANCE 1.0e-20
#endif

//----------------------------------------------------------------------

struct fluxes
{
  long_int LeftFluxStartGlobalIndex [MAX_DIMENSION][MAX_DIMENSION];
  long_int LeftFluxEndGlobalIndex   [MAX_DIMENSION][MAX_DIMENSION];
  long_int RightFluxStartGlobalIndex[MAX_DIMENSION][MAX_DIMENSION];
  long_int RightFluxEndGlobalIndex  [MAX_DIMENSION][MAX_DIMENSION];
  enzo_float *LeftFluxes [MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
  enzo_float *RightFluxes[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION];
};

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

#include "mesh.hpp"
#include "performance.hpp"
#include "method.hpp"
#include "simulation.hpp"
#include "parallel.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "enzo_EnzoNamespace.hpp"

#include "enzo_EnzoFactory.hpp"

#include "enzo_EnzoSimulation.hpp"
#   include "enzo_EnzoSimulationMpi.hpp"
#   include "enzo_EnzoSimulationCharm.hpp"

#include "enzo_EnzoBlock.hpp"

#include "enzo_IoEnzoBlock.hpp"

#include "enzo_EnzoTimestep.hpp"
#include "enzo_EnzoBoundary.hpp"
#include "enzo_EnzoInitialImplosion2.hpp"
#include "enzo_EnzoMethodPpm.hpp"
#include "enzo_EnzoMethodPpml.hpp"

#endif /* ENZO_HPP */

