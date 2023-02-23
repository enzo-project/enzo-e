// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSIdeal.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of EnzoEOSIdeal

#include "cello.hpp"
#include "enzo.hpp"

// #define DEBUG_MATCHING_ARRAY_SHAPES

//----------------------------------------------------------------------

void EnzoEOSIdeal::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change
  PUP::able::pup(p);
  p|gamma_;
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::primitive_from_integration
(const EnzoEFltArrayMap &integration_map, EnzoEFltArrayMap &primitive_map,
 const int stale_depth, const str_vec_t &passive_list,
 const bool ignore_grackle) const
{

  const CelloArray<const enzo_float, 3> density = integration_map.at("density");
  const int mz = density.shape(0);
  const int my = density.shape(1);
  const int mx = density.shape(2);

  // The EOS object doesn't necessarily know what the integration quantities
  // are. This means we take something of an exhaustive approach. This could be
  // more clever if this operations were made a part of the Hydro Integrator
  std::vector<std::string> quantity_list =
    EnzoCenteredFieldRegistry::get_registered_quantities(true, true);

  for (const std::string& key : quantity_list){
    if ((!integration_map.contains(key)) || (!primitive_map.contains(key))){
      continue;
    }

    const CelloArray<const enzo_float, 3> integ_array = integration_map.at(key);
    const CelloArray<enzo_float, 3> prim_array = primitive_map.at(key);

#ifdef DEBUG_MATCHING_ARRAY_SHAPES
    ASSERT6("EnzoEOSIdeal::primitive_from_integration",
            ("The array being copied from integration_map has shape "
             "(%d,%d,%d), while the destination array has shape (%d,%d,%d). "
             "They should be the same."),
            mz, my, mx,
            prim_array.shape(0), prim_array.shape(1), prim_array.shape(2),
            ((prim_array.shape(0) == mz) &&
             (prim_array.shape(1) == my) &&
             (prim_array.shape(2) == mx)));
#endif

    for (int iz = stale_depth; iz < mz - stale_depth; iz++) {
      for (int iy = stale_depth; iy < my - stale_depth; iy++) {
        for (int ix = stale_depth; ix < mx - stale_depth; ix++) {
          prim_array(iz,iy,ix) = integ_array(iz,iy,ix);
        }
      }
    }
  }

  // Convert the passive scalars from conserved-form (i.e. a density) to
  // specific-form (mass fractions)
  for (const std::string& key : passive_list){
    const CelloArray<const enzo_float, 3> cur_conserved
      = integration_map.at(key);
    const CelloArray<enzo_float, 3> out_specific = primitive_map.at(key);

    for (int iz = stale_depth; iz < mz - stale_depth; iz++) {
      for (int iy = stale_depth; iy < my - stale_depth; iy++) {
        for (int ix = stale_depth; ix < mx - stale_depth; ix++) {
          out_specific(iz,iy,ix) = cur_conserved(iz,iy,ix)/density(iz,iy,ix);
        }
      }
    }
  }

  pressure_from_integration(integration_map, primitive_map.at("pressure"),
                            stale_depth, ignore_grackle);
}



//----------------------------------------------------------------------

void EnzoEOSIdeal::pressure_from_integration
(const EnzoEFltArrayMap &integration_map, const EFlt3DArray &pressure,
 const int stale_depth, const bool ignore_grackle) const
{

  const bool mhd = (integration_map.contains("bfield_x") ||
                    integration_map.contains("bfield_y") ||
                    integration_map.contains("bfield_z"));

  const bool idual = enzo::fluid_props()->dual_energy_config().any_enabled();

  EnzoComputePressure::compute_pressure(EnzoFieldAdaptor(integration_map),
                                        pressure, mhd, idual, get_gamma(),
                                        stale_depth, ignore_grackle);
}
