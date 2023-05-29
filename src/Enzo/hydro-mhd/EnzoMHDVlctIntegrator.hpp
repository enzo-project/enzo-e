// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMHDVlctIntegrator.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 28 2023
/// @brief    [\ref Enzo] Declaration of EnzoMHDVlctIntegrator

#ifndef ENZO_ENZO_MHD_VLCT_INTEGRATOR_HPP
#define ENZO_ENZO_MHD_VLCT_INTEGRATOR_HPP

class EnzoMHDVlctIntegrator {

  /// @class    EnzoMHDVlctIntegrator
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Solve equations for VL + CT MHD method

public:

  EnzoMHDVlctIntegrator
  (EnzoRiemann* riemann_solver,
   EnzoReconstructor* half_dt_recon,
   EnzoReconstructor* full_dt_recon,
   EnzoIntegrationQuanUpdate *integration_quan_updater)
    : riemann_solver_(riemann_solver),
      half_dt_recon_(half_dt_recon),
      full_dt_recon_(full_dt_recon),
      integration_quan_updater_(integration_quan_updater)
  { }

  ~EnzoMHDVlctIntegrator();

public:
  /// Pointer to the Riemann solver
  EnzoRiemann *riemann_solver_;
  /// Pointer to the reconstructor used to reconstruct the fluid during the
  /// first half time-step (usually nearest-neighbor)
  EnzoReconstructor *half_dt_recon_;
  /// Pointer to the reconstructor used to reconstruct the fluid during the
  /// full time-step
  EnzoReconstructor *full_dt_recon_;
  /// Pointer to the integration quantity updater
  EnzoIntegrationQuanUpdate *integration_quan_updater_;

};





#endif /* ENZO_ENZO_MHD_VLCT_INTEGRATOR_HPP */
