// See LICENSE_CELLO file for license and copyright information

/// @file     enzo.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Include file for the \ref Enzo component

#ifndef ENZO_PRIVATE_HPP
#define ENZO_PRIVATE_HPP

//----------------------------------------------------------------------

#define OMEGA_TOLERANCE 1.0e-5
 
#ifdef CONFIG_PRECISION_SINGLE
#   define ETA_TOLERANCE 1.0e-5
#   define ENZO_HUGE_VAL HUGE_VALF
#   define scalar_type_enzo_float scalar_type_float
#   define CONFIG_BFLOAT_4
#   define CONFIG_PFLOAT_4
#endif
#ifdef CONFIG_PRECISION_DOUBLE
#   define ETA_TOLERANCE 1.0e-10
#   define ENZO_HUGE_VAL HUGE_VAL
#   define scalar_type_enzo_float scalar_type_double
#   define CONFIG_BFLOAT_8
#   define CONFIG_PFLOAT_8
#endif
#ifdef CONFIG_PRECISION_QUAD
#   define ETA_TOLERANCE 1.0e-20
#   define ENZO_HUGE_VAL HUGE_VALL
#   define CONFIG_BFLOAT_16
#   define CONFIG_PFLOAT_16
#endif

//----------------------------------------------------------------------

#define INDEX_TURBULENCE_VAD   0
#define INDEX_TURBULENCE_AAD   1
#define INDEX_TURBULENCE_VVDoT 2
#define INDEX_TURBULENCE_VVoT  3
#define INDEX_TURBULENCE_VVD   4
#define INDEX_TURBULENCE_VV    5
#define INDEX_TURBULENCE_DD    6

#define INDEX_TURBULENCE_DAx   7
#define INDEX_TURBULENCE_DAy   8
#define INDEX_TURBULENCE_DAz   9

#define INDEX_TURBULENCE_DVx   10
#define INDEX_TURBULENCE_DVy   11
#define INDEX_TURBULENCE_DVz   12

#define INDEX_TURBULENCE_DlnD  13

/* minD and maxD must be the last two */
#define INDEX_TURBULENCE_minD  14
#define INDEX_TURBULENCE_maxD  15

#define MAX_TURBULENCE_ARRAY 16 /* size of global array of global values to track */

//----------------------------------------------------------------------

// #include "macros_and_parameters.h"
#include "enzo_defines.hpp"
#include "enzo_typedefs.hpp"
#include "enzo_fortran.hpp"

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

enum hydro_type 
  {
    hydro_unknown,
    hydro_ppm,
    hydro_ppml
  };

//----------------------------------------------------------------------

enum return_enum {
  return_unknown,
  return_converged,
  return_error_max_iter_reached
};

//----------------------------------------------------------------------
// WARNING 100 must be larger than number of phases in
// src/Cello/_mesh.hpp phase_type

enum enzo_phase_type 
  {
    phase_enzo_first = 100,
    phase_enzo_matvec = phase_enzo_first
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

#ifdef CONFIG_USE_GRACKLE
#   include "grackle.h"
#endif

//----------------------------------------------------------------------

#include "fortran.h" /* included so scons knowns to install fortran.h */

#include "enzo_EnzoFactory.hpp"

#include "enzo_EnzoSimulation.hpp"

#include "enzo_EnzoProblem.hpp"

#include "enzo_EnzoConfig.hpp"

#include "enzo_EnzoBlock.hpp"

#include "enzo_IoEnzoBlock.hpp"

#include "enzo_EnzoBoundary.hpp"

#include "enzo_EnzoInitialGrackleTest.hpp"
#include "enzo_EnzoInitialImplosion2.hpp"
#include "enzo_EnzoInitialSedovArray2.hpp"
#include "enzo_EnzoInitialSedovArray3.hpp"
#include "enzo_EnzoInitialTurbulence.hpp"

#include "enzo_EnzoRefineShock.hpp"

#include "enzo_EnzoMethodNull.hpp"
#include "enzo_EnzoMethodPpm.hpp"
#include "enzo_EnzoMethodPpml.hpp"
#include "enzo_EnzoMethodHeat.hpp"
#include "enzo_EnzoMethodGrackle.hpp"
#include "enzo_EnzoMethodTurbulence.hpp"
#include "enzo_EnzoMethodGravityCg.hpp"
#include "enzo_EnzoMethodGravityBiCGStab.hpp"

#include "enzo_EnzoMatrixLaplace.hpp"
#include "enzo_EnzoMatrixDiagonal.hpp"
#include "enzo_EnzoMatrixIdentity.hpp"

#include "enzo_EnzoComputePressure.hpp"
#include "enzo_EnzoComputeTemperature.hpp"
#include "enzo_EnzoComputeAcceleration.hpp"

#include "enzo_EnzoProlong.hpp"
#include "enzo_EnzoProlongMC1.hpp"
#include "enzo_EnzoProlongPoisson.hpp"
#include "enzo_EnzoRestrict.hpp"

#endif /* ENZO_PRIVATE_HPP */

