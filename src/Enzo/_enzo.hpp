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
// Cello include file
//----------------------------------------------------------------------

#include "cello.hpp"

//----------------------------------------------------------------------
// Component class includes
//----------------------------------------------------------------------

#include "fortran.h" /* included so scons knowns to install fortran.h */

#include "fortran_types.h" /* included so scons knowns to install fortran.h */

#include "enzo_constants.hpp"

// TODO: remove this include statement (this may be easier to do once we have
//       finished separating out the various subcomponents)
#include "utils/utils.hpp"

#include "inference/Index3.hpp"

// TODO: remove the following 3 after factoring out other subcomponents
#include "Enzo/cosmology/cosmology.hpp"
#include "Enzo/chemistry/chemistry.hpp"
#include "Enzo/fluid-props/fluid-props.hpp"

#include "enzo-core/EnzoUnits.hpp"

#include "enzo-core/EnzoFactory.hpp"

#include "enzo-core/EnzoSimulation.hpp"

#include "enzo-core/EnzoProblem.hpp"

#include "enzo-core/EnzoConfig.hpp"
#include "enzo-core/EnzoState.hpp"
#include "enzo-core/EnzoBlock.hpp"

#include "inference/EnzoLevelArray.hpp"
#include "inference/EnzoMethodInference.hpp"

#include "enzo-core/EnzoBoundary.hpp"

#include "enzo-core/EnzoMethodBalance.hpp"

#include "enzo-core/EnzoMsgCheck.hpp"

#include "enzo-core/EnzoStopping.hpp"

#endif /* ENZO_PRIVATE_HPP */
