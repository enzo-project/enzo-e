// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSourceInternalEnergy.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 6 2022
/// @brief    [\ref Enzo] Implementation of Enzo's SourceGravity class
///
/// This is adapted from a snippet of hydro_rk/Grid_SourceTerms.C in the
/// original Enzo codebase

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoSourceGravity::calculate_source
(const double dt, const EnzoEFltArrayMap &prim_map,
 EnzoEFltArrayMap &dUcons_map, const EnzoEFltArrayMap &accel_map,
 const int stale_depth) const noexcept
{
  // SANITY CHECK:
  ASSERT("EnzoSourceInternalEnergy::calculate_source",
	 "Barotropic equations of state are not currently supported.",
	 !(enzo::fluid_props()->has_barotropic_eos()) );

  // load primitives
  const CelloView<const enzo_float,3> density = prim_map.at("density");
  const CelloView<const enzo_float,3> velocity_x = prim_map.at("velocity_x");
  const CelloView<const enzo_float,3> velocity_y = prim_map.at("velocity_y");
  const CelloView<const enzo_float,3> velocity_z = prim_map.at("velocity_z");

  // load acceleration fields
  const CelloView<const enzo_float,3> accel_x = accel_map.at("acceleration_x");
  const CelloView<const enzo_float,3> accel_y = accel_map.at("acceleration_y");
  const CelloView<const enzo_float,3> accel_z = accel_map.at("acceleration_z");

  // load accumulation arrays (reminder: these correspond to the conserved
  // forms of the integration quantities)
  const CelloView<enzo_float,3> dmom_x = dUcons_map.at("velocity_x");
  const CelloView<enzo_float,3> dmom_y = dUcons_map.at("velocity_y");
  const CelloView<enzo_float,3> dmom_z = dUcons_map.at("velocity_z");
  const CelloView<enzo_float,3> dE_dens = dUcons_map.at("total_energy");

  const int mz = density.shape(0);
  const int my = density.shape(1);
  const int mx = density.shape(2);

  for (int iz = stale_depth; iz < (mz - stale_depth); iz++){
    for (int iy = stale_depth; iy < (my - stale_depth); iy++){
      for (int ix = stale_depth; ix < (mx - stale_depth); ix++){

        enzo_float rho = density(iz,iy,ix);
        enzo_float ax = accel_x(iz,iy,ix);
        enzo_float ay = accel_y(iz,iy,ix);
        enzo_float az = accel_z(iz,iy,ix);

        dmom_x(iz,iy,ix) += dt * rho * ax;
        dmom_y(iz,iy,ix) += dt * rho * ay;
        dmom_z(iz,iy,ix) += dt * rho * az;
        dE_dens(iz,iy,ix) += dt * rho * ( (velocity_x(iz,iy,ix) * ax) +
                                          (velocity_y(iz,iy,ix) * ay) +
                                          (velocity_z(iz,iy,ix) * az));

      }
    }
  }

}
