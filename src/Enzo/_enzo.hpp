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
  enzo_sync_id_method_grackle,
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
  #define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC
  #include <grackle.h>
}
#else
extern "C" { // declare the names of Grackle types so can reduce the usage of
             // ifdef statements
  struct chemistry_data;
  struct chemistry_data_storage;
  struct code_units;
  struct grackle_field_data;
}
#endif

//----------------------------------------------------------------------

#include "fortran.h" /* included so scons knowns to install fortran.h */

#include "fortran_types.h" /* included so scons knowns to install fortran.h */

#include "enzo_constants.hpp"
#include "utils/EnzoEFltArrayMap.hpp"
#include "utils/utils.hpp"

#include "enzo_Index3.hpp"

#include "cosmology/EnzoPhysicsCosmology.hpp"

// [order dependencies:]
#include "fluid-props/EnzoEOSIdeal.hpp"
#include "fluid-props/EnzoEOSIsothermal.hpp"
#include "fluid-props/EnzoEOSVariant.hpp"

#include "fluid-props/EnzoDualEnergyConfig.hpp"
#include "fluid-props/EnzoFluidFloorConfig.hpp"
#include "fluid-props/EnzoPhysicsFluidProps.hpp"

#include "enzo-core/EnzoUnits.hpp"

// [order dependencies:]
#include "utils/EnzoEFltArrayMap.hpp"
#include "utils/EnzoFieldAdaptor.hpp"
#include "chemistry/GrackleChemistryData.hpp"
#include "chemistry/GrackleFacade.hpp"

#include "enzo-core/EnzoFactory.hpp"

#include "enzo-core/EnzoSimulation.hpp"

#include "enzo-core/EnzoProblem.hpp"

#include "enzo-core/EnzoConfig.hpp"

#include "enzo-core/EnzoBlock.hpp"


#include "enzo_EnzoLevelArray.hpp"

#include "io/IoEnzoBlock.hpp"
#include "io/IoEnzoReader.hpp"
#include "io/IoEnzoWriter.hpp"

#include "enzo-core/EnzoBoundary.hpp"

#include "initial/EnzoInitialBCenter.hpp"
#include "initial/EnzoInitialCloud.hpp"
#include "initial/EnzoInitialCollapse.hpp"
#include "initial/EnzoInitialCosmology.hpp"
#include "tests/EnzoInitialFeedbackTest.hpp"
#include "tests/EnzoInitialGrackleTest.hpp"
#include "io/EnzoInitialHdf5.hpp"
#include "initial/EnzoInitialImplosion2.hpp"
#include "initial/EnzoInitialInclinedWave.hpp"
#include "io/EnzoInitialMusic.hpp"
#include "obsolete/EnzoInitialPm.hpp"
#include "tests/EnzoInitialPpmlTest.hpp"
#include "initial/EnzoInitialSedovArray2.hpp"
#include "initial/EnzoInitialSedovArray3.hpp"
#include "initial/EnzoInitialSedovRandom.hpp"
#include "initial/EnzoInitialShockTube.hpp"
#include "initial/EnzoInitialSoup.hpp"
#include "initial/EnzoInitialTurbulence.hpp"
#include "initial/EnzoInitialIsolatedGalaxy.hpp"
#include "initial/EnzoInitialBurkertBodenheimer.hpp"
#include "tests/EnzoInitialMergeSinksTest.hpp"
#include "tests/EnzoInitialAccretionTest.hpp"
#include "initial/EnzoInitialShuCollapse.hpp"
#include "tests/EnzoInitialBBTest.hpp"

#include "mesh/EnzoRefineShock.hpp"
#include "mesh/EnzoRefineParticleMass.hpp"
#include "mesh/EnzoRefineMass.hpp"

// [order dependencies:]
#include "particle/formation/EnzoSinkParticle.hpp"
#include "particle/formation/EnzoBondiHoyleSinkParticle.hpp"
#include "particle/formation/EnzoFluxSinkParticle.hpp"

#include "utils/EnzoCenteredFieldRegistry.hpp"
#include "hydro-mhd/EnzoIntegrationQuanUpdate.hpp"
#include "hydro-mhd/EnzoLazyPassiveScalarFieldList.hpp"
#include "hydro-mhd/EnzoPermutedCoordinates.hpp"
#include "hydro-mhd/EnzoSourceGravity.hpp"
#include "hydro-mhd/EnzoSourceInternalEnergy.hpp"

#include "hydro-mhd/EnzoReconstructor.hpp"
#include "hydro-mhd/EnzoReconstructorNN.hpp"
#include "hydro-mhd/EnzoReconstructorPLM.hpp"

// public header for the EnzoRiemann sub-library. This needs to be included
// after the headers for:
//     EnzoEFltArrayMap, EnzoCenteredFieldRegistry, & EnzoEquationOfState
// but before the header for EnzoMethodMHDVlct EnzoBfieldMethod and EnzoBfieldMethodCT
#include "hydro-mhd/riemann/EnzoRiemann.hpp"

// [order dependencies:]
#include "hydro-mhd/EnzoBfieldMethod.hpp"
#include "hydro-mhd/EnzoBfieldMethodCT.hpp"


#include "enzo_EnzoMethodInference.hpp"

#include "particle/formation/EnzoMethodAccretion.hpp"
#include "gravity/EnzoMethodBackgroundAcceleration.hpp"
#include "particle/formation/EnzoMethodBondiHoyleAccretion.hpp"
#include "enzo-core/EnzoMethodBalance.hpp"
#include "io/EnzoMethodCheck.hpp"
#include "cosmology/EnzoMethodComovingExpansion.hpp"
#include "cosmology/EnzoMethodCosmology.hpp"
#include "particle/feedback/EnzoMethodDistributedFeedback.hpp"
#include "particle/feedback/EnzoMethodFeedback.hpp"
#include "particle/feedback/EnzoMethodFeedbackSTARSS.hpp"
#include "particle/formation/EnzoMethodFluxAccretion.hpp"
#include "chemistry/EnzoMethodGrackle.hpp"
#include "gravity/EnzoMethodGravity.hpp"
#include "assorted/EnzoMethodHeat.hpp"
#include "obsolete/EnzoMethodHydro.hpp"
#include "particle/formation/EnzoMethodMergeSinks.hpp"
#include "hydro-mhd/EnzoMethodMHDVlct.hpp"
#include "assorted/EnzoMethodM1Closure.hpp"
#include "gravity/EnzoMethodPmDeposit.hpp"
#include "particle/EnzoMethodPmUpdate.hpp"
#include "hydro-mhd/EnzoMethodPpm.hpp"
#include "hydro-mhd/EnzoMethodPpml.hpp"
#include "particle/formation/EnzoMethodSinkMaker.hpp"
#include "particle/formation/EnzoMethodStarMaker.hpp"
#include "particle/formation/EnzoMethodStarMakerSTARSS.hpp"
#include "particle/formation/EnzoMethodStarMakerStochasticSF.hpp"
#include "particle/formation/EnzoMethodThresholdAccretion.hpp"
#include "assorted/EnzoMethodTurbulence.hpp"

#include "gravity/matrix/EnzoMatrixDiagonal.hpp"
#include "gravity/matrix/EnzoMatrixIdentity.hpp"
#include "gravity/matrix/EnzoMatrixLaplace.hpp"

#include "enzo-core/EnzoMsgCheck.hpp"

#include "gravity/EnzoComputeAcceleration.hpp"
#include "utils/EnzoComputeCicInterp.hpp"
#include "fluid-props/EnzoComputePressure.hpp"
#include "fluid-props/EnzoComputeTemperature.hpp"
#include "chemistry/EnzoComputeCoolingTime.hpp"

#include "gravity/solvers/EnzoSolverBiCgStab.hpp"
#include "gravity/solvers/EnzoSolverCg.hpp"
#include "gravity/solvers/EnzoSolverDd.hpp"
#include "gravity/solvers/EnzoSolverDiagonal.hpp"
#include "gravity/solvers/EnzoSolverJacobi.hpp"
#include "gravity/solvers/EnzoSolverMg0.hpp"

#include "enzo-core/EnzoStopping.hpp"

#include "mesh/EnzoProlong.hpp"
#include "mesh/EnzoProlongMC1.hpp"
#include "mesh/EnzoProlongPoisson.hpp"
#include "mesh/EnzoRestrict.hpp"

#endif /* ENZO_PRIVATE_HPP */
