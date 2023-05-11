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
#   define type_enzo_float type_float
#endif
#ifdef CONFIG_PRECISION_DOUBLE
#   define ETA_TOLERANCE 1.0e-10
#   define ENZO_HUGE_VAL HUGE_VAL
#   define type_enzo_float type_double
#endif
#ifdef CONFIG_PRECISION_QUAD
#   define ETA_TOLERANCE 1.0e-20
#   define ENZO_HUGE_VAL HUGE_VALL
#endif

//----------------------------------------------------------------------

enum mass_type {
  mass_unknown,
  mass_dark,
  mass_baryon
};

//----------------------------------------------------------------------

enum {
  index_turbulence_vad,
  index_turbulence_aad,
  index_turbulence_vvdot,
  index_turbulence_vvot,
  index_turbulence_vvd,
  index_turbulence_vv,
  index_turbulence_dd,
  index_turbulence_d,
  index_turbulence_dax,
  index_turbulence_day,
  index_turbulence_daz,
  index_turbulence_dvx,
  index_turbulence_dvy,
  index_turbulence_dvz,
  index_turbulence_dlnd,
  index_turbulence_zones,
  index_turbulence_mind,
  index_turbulence_maxd,
  max_turbulence_array };

//----------------------------------------------------------------------

enum enzo_sync_id {
  enzo_sync_id_cg = sync_id_last,
  enzo_sync_id_comoving_expansion,
  enzo_sync_id_method_background_acceleration,
  enzo_sync_id_method_cosmology,
  enzo_sync_id_method_feedback,
  enzo_sync_id_method_radiative_transfer,
#ifdef CONFIG_USE_GRACKLE
  enzo_sync_id_method_grackle,
#endif
  enzo_sync_id_method_gravity,
  enzo_sync_id_method_gravity_continue,
  enzo_sync_id_method_heat,
  enzo_sync_id_method_null,
  enzo_sync_id_method_pm_deposit,
  enzo_sync_id_method_pm_update,
  enzo_sync_id_method_ppm,
  enzo_sync_id_method_ppml,
  enzo_sync_id_method_star_maker,
  enzo_sync_id_method_turbulence,
  enzo_sync_id_method_vlct,
  enzo_sync_id_solver_bicgstab,
  enzo_sync_id_solver_bicgstab_precon_1,
  enzo_sync_id_solver_bicgstab_precon_2,
  enzo_sync_id_solver_bicgstab_loop_25,
  enzo_sync_id_solver_bicgstab_loop_85,
  enzo_sync_id_solver_cg_matvec,
  enzo_sync_id_solver_cg_loop_2,
  enzo_sync_id_solver_dd,
  enzo_sync_id_solver_dd_coarse,
  enzo_sync_id_solver_dd_domain,
  enzo_sync_id_solver_dd_smooth,
  enzo_sync_id_solver_mg0,
  enzo_sync_id_solver_mg0_coarse,
  enzo_sync_id_solver_mg0_last,
  enzo_sync_id_solver_mg0_post,
  enzo_sync_id_solver_mg0_pre,
  enzo_sync_id_solver_jacobi_1,
  enzo_sync_id_solver_jacobi_2,
  enzo_sync_id_solver_jacobi_3
};

//----------------------------------------------------------------------

// #include "macros_and_parameters.h"
#include "enzo_defines.hpp"
#include "enzo_typedefs.hpp"
#include "enzo_fortran.hpp"
#include "enzo_reductions.hpp"

//----------------------------------------------------------------------

enum return_enum {
  return_unknown,
  return_converged,
  return_diverged,
  return_error,
  return_bypass
};

//----------------------------------------------------------------------

const int field_undefined = -1;

//----------------------------------------------------------------------

struct enzo_fluxes
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
#include <stdlib.h>
extern "C" {
  #include <grackle.h>
}
#endif

//----------------------------------------------------------------------

#include "fortran.h" /* included so scons knowns to install fortran.h */

#include "fortran_types.h" /* included so scons knowns to install fortran.h */

#include "enzo_constants.hpp"

#include "enzo_EnzoPhysicsCosmology.hpp"

#include "enzo_EnzoDualEnergyConfig.hpp"
#include "enzo_EnzoFluidFloorConfig.hpp"
#include "enzo_EnzoPhysicsFluidProps.hpp"

#include "enzo_EnzoUnits.hpp"

#include "enzo_EnzoFactory.hpp"

#include "enzo_EnzoSimulation.hpp"

#include "enzo_EnzoProblem.hpp"

#include "enzo_EnzoConfig.hpp"

#include "enzo_EnzoBlock.hpp"

#include "enzo_IoEnzoBlock.hpp"
#include "enzo_IoEnzoReader.hpp"
#include "enzo_IoEnzoWriter.hpp"

#include "enzo_EnzoBoundary.hpp"

#include "enzo_EnzoInitialBCenter.hpp"
#include "enzo_EnzoInitialCloud.hpp"
#include "enzo_EnzoInitialCollapse.hpp"
#include "enzo_EnzoInitialCosmology.hpp"
#include "enzo_EnzoInitialFeedbackTest.hpp"
#include "enzo_EnzoInitialGrackleTest.hpp"
#include "enzo_EnzoInitialHdf5.hpp"
#include "enzo_EnzoInitialImplosion2.hpp"
#include "enzo_EnzoInitialInclinedWave.hpp"
#include "enzo_EnzoInitialMusic.hpp"
#include "enzo_EnzoInitialPm.hpp"
#include "enzo_EnzoInitialPpmlTest.hpp"
#include "enzo_EnzoInitialSedovArray2.hpp"
#include "enzo_EnzoInitialSedovArray3.hpp"
#include "enzo_EnzoInitialSedovRandom.hpp"
#include "enzo_EnzoInitialShockTube.hpp"
#include "enzo_EnzoInitialSoup.hpp"
#include "enzo_EnzoInitialTurbulence.hpp"
#include "enzo_EnzoInitialIsolatedGalaxy.hpp"
#include "enzo_EnzoInitialBurkertBodenheimer.hpp"
#include "enzo_EnzoInitialMergeSinksTest.hpp"
#include "enzo_EnzoInitialAccretionTest.hpp"
#include "enzo_EnzoInitialShuCollapse.hpp"
#include "enzo_EnzoInitialBBTest.hpp"

#include "enzo_EnzoRefineShock.hpp"
#include "enzo_EnzoRefineParticleMass.hpp"
#include "enzo_EnzoRefineMass.hpp"

// [order dependencies:]
#include "enzo_EnzoSinkParticle.hpp"
#include "enzo_EnzoBondiHoyleSinkParticle.hpp"
#include "enzo_EnzoFluxSinkParticle.hpp"

// [order dependencies:]
#include "enzo_EnzoEFltArrayMap.hpp"
#include "enzo_EnzoEquationOfState.hpp"
#include "enzo_EnzoEOSIdeal.hpp"

#include "enzo_EnzoCenteredFieldRegistry.hpp"
#include "enzo_EnzoFieldAdaptor.hpp"
#include "enzo_EnzoIntegrationQuanUpdate.hpp"
#include "enzo_EnzoLazyPassiveScalarFieldList.hpp"
#include "enzo_EnzoPermutedCoordinates.hpp"
#include "enzo_EnzoReconstructor.hpp"
#include "enzo_EnzoReconstructorNN.hpp"
#include "enzo_EnzoReconstructorPLM.hpp"
#include "enzo_EnzoSourceGravity.hpp"
#include "enzo_EnzoSourceInternalEnergy.hpp"

// public header for the EnzoRiemann sub-library. This needs to be included
// after the headers for:
//     EnzoEFltArrayMap, EnzoCenteredFieldRegistry, & EnzoEquationOfState
// but before the header for EnzoMethodMHDVlct EnzoBfieldMethod and EnzoBfieldMethodCT
#include "EnzoRiemann.hpp"

// [order dependencies:]
#include "enzo_EnzoBfieldMethod.hpp"
#include "enzo_EnzoBfieldMethodCT.hpp"

#include "enzo_EnzoMethodAccretion.hpp"
#include "enzo_EnzoMethodBackgroundAcceleration.hpp"
#include "enzo_EnzoMethodBondiHoyleAccretion.hpp"
#include "enzo_EnzoMethodCheck.hpp"
#include "enzo_EnzoMethodComovingExpansion.hpp"
#include "enzo_EnzoMethodCosmology.hpp"
#include "enzo_EnzoMethodDistributedFeedback.hpp"
#include "enzo_EnzoMethodFeedback.hpp"
#include "enzo_EnzoMethodFeedbackSTARSS.hpp"
#include "enzo_EnzoMethodFluxAccretion.hpp"
#include "enzo_EnzoMethodGrackle.hpp"
#include "enzo_EnzoMethodGravity.hpp"
#include "enzo_EnzoMethodHeat.hpp"
#include "enzo_EnzoMethodHydro.hpp"
#include "enzo_EnzoMethodMergeSinks.hpp"
#include "enzo_EnzoMethodMHDVlct.hpp"
#include "enzo_EnzoMethodPmDeposit.hpp"
#include "enzo_EnzoMethodPmUpdate.hpp"
#include "enzo_EnzoMethodPpm.hpp"
#include "enzo_EnzoMethodPpml.hpp"
#include "enzo_EnzoMethodBalance.hpp"
#include "enzo_EnzoMethodSinkMaker.hpp"
#include "enzo_EnzoMethodStarMaker.hpp"
#include "enzo_EnzoMethodStarMakerSTARSS.hpp"
#include "enzo_EnzoMethodStarMakerStochasticSF.hpp"
#include "enzo_EnzoMethodMHDVlct.hpp"
#include "enzo_EnzoMethodM1Closure.hpp"
#include "enzo_EnzoMethodThresholdAccretion.hpp"
#include "enzo_EnzoMethodTurbulence.hpp"

#include "enzo_EnzoMatrixDiagonal.hpp"
#include "enzo_EnzoMatrixIdentity.hpp"
#include "enzo_EnzoMatrixLaplace.hpp"

#include "enzo_EnzoMsgCheck.hpp"

#include "enzo_EnzoComputeAcceleration.hpp"
#include "enzo_EnzoComputeCicInterp.hpp"
#include "enzo_EnzoComputePressure.hpp"
#include "enzo_EnzoComputeTemperature.hpp"
#ifdef CONFIG_USE_GRACKLE
  #include "enzo_EnzoComputeCoolingTime.hpp"
#endif

#include "enzo_EnzoSolverBiCgStab.hpp"
#include "enzo_EnzoSolverCg.hpp"
#include "enzo_EnzoSolverDd.hpp"
#include "enzo_EnzoSolverDiagonal.hpp"
#include "enzo_EnzoSolverJacobi.hpp"
#include "enzo_EnzoSolverMg0.hpp"

#include "enzo_EnzoStopping.hpp"

#include "enzo_EnzoProlong.hpp"
#include "enzo_EnzoProlongMC1.hpp"
#include "enzo_EnzoProlongPoisson.hpp"
#include "enzo_EnzoRestrict.hpp"

#endif /* ENZO_PRIVATE_HPP */
