// See LICENSE_CELLO file for license and copyright information

/// @file     enzo.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Include file for the \ref Enzo component

#ifndef ENZO_PRIVATE_HPP
#define ENZO_PRIVATE_HPP

//----------------------------------------------------------------------

#include "enzo_defines.hpp"
#include "enzo_typedefs.hpp"
#include "enzo_fortran.hpp"

//----------------------------------------------------------------------

#define OMEGA_TOLERANCE 1.0e-5
 
#ifdef CONFIG_PRECISION_SINGLE
#   define ETA_TOLERANCE 1.0e-5
#   define ENZO_HUGE_VAL HUGE_VALF
#   define scalar_type_enzo_float scalar_type_float
#endif
#ifdef CONFIG_PRECISION_DOUBLE
#   define ETA_TOLERANCE 1.0e-10
#   define ENZO_HUGE_VAL HUGE_VAL
#   define scalar_type_enzo_float scalar_type_double
#endif
#ifdef CONFIG_PRECISION_QUAD
#   define ETA_TOLERANCE 1.0e-20
#   define ENZO_HUGE_VAL HUGE_VALL
#endif

//----------------------------------------------------------------------

enum bc_enum 
  { // explicitly enumerated to match what Enzo expects
    bc_unknown    = 0, 
    bc_reflecting = 1, 
    bc_outflow    = 2, 
    bc_inflow     = 3, 
    bc_periodic   = 4 
  };

//----------------------------------------------------------------------

enum hydro_type {
  hydro_unknown,
  hydro_ppm,
  hydro_ppml
};

//----------------------------------------------------------------------

  const int field_undefined = -1;

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
// Cello include file
//----------------------------------------------------------------------

#include "cello.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "grackle.h"

//----------------------------------------------------------------------

#include "fortran.h" /* included so scons knowns to install fortran.h */

#include "enzo_EnzoFactory.hpp"

#include "enzo_EnzoSimulationMpi.hpp"
#include "enzo_EnzoSimulationCharm.hpp"

#include "enzo_EnzoProblem.hpp"

#include "enzo_EnzoConfig.hpp"

#include "enzo_EnzoBlock.hpp"

#include "enzo_IoEnzoBlock.hpp"

// #include "enzo_EnzoTimestep.hpp"
// #include "enzo_EnzoTimestepPpml.hpp"
#include "enzo_EnzoBoundary.hpp"

#include "enzo_EnzoInitialImplosion2.hpp"
#include "enzo_EnzoInitialSedovArray2.hpp"
#include "enzo_EnzoInitialSedovArray3.hpp"

#include "enzo_EnzoMethodNull.hpp"
#include "enzo_EnzoMethodPpm.hpp"
#include "enzo_EnzoMethodPpml.hpp"
#include "enzo_EnzoMethodHeat.hpp"
#include "enzo_EnzoMethodGrackle.hpp"

#include "enzo_EnzoProlong.hpp"
#include "enzo_EnzoProlongMC1.hpp"
#include "enzo_EnzoRestrict.hpp"


#endif /* ENZO_PRIVATE_HPP */

