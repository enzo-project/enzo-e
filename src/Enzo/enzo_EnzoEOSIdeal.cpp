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
  p|density_floor_;
  p|pressure_floor_;
  p|dual_energy_formalism_;
  p|dual_energy_formalism_eta_;
}

//----------------------------------------------------------------------

// returns true if grackle is in use and if gamma can vary spatially
bool grackle_variable_gamma_(){
#ifdef CONFIG_USE_GRACKLE
  if (enzo::config()->method_grackle_use_grackle){
    if (enzo::config()->method_grackle_chemistry->primordial_chemistry > 1) {
      return true;
    }
  }
#endif
  return false;
}

//----------------------------------------------------------------------

// Helper function that performs a quick check to confirm that certain fields
// are held by both reconstructable_group and integrable_group
void confirm_same_kv_pair_(const EnzoEFltArrayMap &reconstructable,
                           const EnzoEFltArrayMap &integrable,
                           bool allowed_to_be_omitted, const std::string &key,
                           const std::string &func_name)
{
  if (allowed_to_be_omitted && ( (!integrable.contains(key)) ||
                                 (!reconstructable.contains(key)) ) ){
    return;
  } else {
    ASSERT1(func_name.c_str(),
            ("The arrays associated with the \"%s\" key in the recontructable "
             "and integrable maps are expected to be aliases of each other."),
            key.c_str(), integrable.at(key).is_alias(reconstructable.at(key)));
  }
}

//----------------------------------------------------------------------

void check_recon_integ_overlap_
(const EnzoEFltArrayMap &reconstructable, const EnzoEFltArrayMap &integrable,
 const std::string &func_name, const str_vec_t &passive_list)
{
  // We assume that the following groups are represented by the same fields in
  // integrable and reconstructable. This should probably not be hardcoded
  const std::array<std::string,4> standard_common_keys =
    {"density", "velocity_x", "velocity_y", "velocity_z"};
  const std::array<std::string,3> optional_common_keys =
    {"bfield_x", "bfield_y", "bfield_z"};

  for (const std::string& key : standard_common_keys){
    confirm_same_kv_pair_(reconstructable, integrable, false, key, func_name);
  }
  for (const std::string& key : optional_common_keys){
    confirm_same_kv_pair_(reconstructable, integrable, true, key, func_name);
  }
  for (const std::string& key : passive_list){
    confirm_same_kv_pair_(reconstructable, integrable, false, key, func_name);
  }
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::reconstructable_from_integrable
(EnzoEFltArrayMap &integrable, EnzoEFltArrayMap &reconstructable,
 EnzoEFltArrayMap &conserved_passive_map, int stale_depth,
 const str_vec_t &passive_list) const
{

  // Confirm that the expected fields (e.g. density, vx, vy, vz, bx, by, bz)
  // are the same in reconstructable_group and integrable_group
  check_recon_integ_overlap_(reconstructable, integrable,
			     "EnzoEOSIdeal::reconstructable_from_integrable",
                             passive_list);

  // Simply compute the pressure
  pressure_from_integrable(integrable, reconstructable.at("pressure"),
                           conserved_passive_map, stale_depth);
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::integrable_from_reconstructable
(EnzoEFltArrayMap &reconstructable, EnzoEFltArrayMap &integrable,
 int stale_depth, const str_vec_t &passive_list) const
{
  if (grackle_variable_gamma_()){
    // we would need to model a field with the
    ERROR("EnzoEOSIdeal::pressure_from_integrable",
	  "Doesn't currently support spatial variations in gamma");
  }

  const bool idual = this->uses_dual_energy_formalism();
  const bool mag   = (reconstructable.contains("bfield_x") ||
                      reconstructable.contains("bfield_y") ||
                      reconstructable.contains("bfield_z"));
  // Confirm that the expected fields (e.g. density, vx, vy, vz, bx, by, bz)
  // are the same in reconstructable_group and integrable_group
  check_recon_integ_overlap_(reconstructable, integrable,
			     "EnzoEOSIdeal::integrable_from_reconstructable",
                             passive_list);

  EFlt3DArray density = reconstructable.get("density", stale_depth);
  EFlt3DArray vx = reconstructable.get("velocity_x", stale_depth);
  EFlt3DArray vy = reconstructable.get("velocity_y", stale_depth);
  EFlt3DArray vz = reconstructable.get("velocity_z", stale_depth);
  EFlt3DArray pressure = reconstructable.get("pressure", stale_depth);

  EFlt3DArray bx, by, bz;
  if (mag){
    bx = reconstructable.get("bfield_x", stale_depth);
    by = reconstructable.get("bfield_y", stale_depth);
    bz = reconstructable.get("bfield_z", stale_depth);
  }

  EFlt3DArray eint, etot;
  if (idual){
    eint = integrable.get("internal_energy", stale_depth);
  }
  etot = integrable.get("total_energy", stale_depth);

  enzo_float inv_gm1 = 1./(get_gamma()-1.);

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	enzo_float v2 = (vx(iz,iy,ix) * vx(iz,iy,ix) +
			 vy(iz,iy,ix) * vy(iz,iy,ix) +
			 vz(iz,iy,ix) * vz(iz,iy,ix));
	enzo_float inv_rho = 1./density(iz,iy,ix);
	enzo_float eint_val = pressure(iz,iy,ix) * inv_gm1 * inv_rho;
	if (idual){
	  eint(iz,iy,ix) = eint_val;
	}
        enzo_float etot_val = eint_val + (0.5 * v2);
        if (mag){
          enzo_float b2 = (bx(iz,iy,ix) * bx(iz,iy,ix) +
                           by(iz,iy,ix) * by(iz,iy,ix) +
                           bz(iz,iy,ix) * bz(iz,iy,ix));
          etot_val += (0.5 * b2 * inv_rho);
        }
        etot(iz,iy,ix) = etot_val;
      }
    }
  }
}



//----------------------------------------------------------------------

void EnzoEOSIdeal::pressure_from_integrable
(EnzoEFltArrayMap &integrable_map, const EFlt3DArray &pressure,
 EnzoEFltArrayMap &conserved_passive_map, int stale_depth) const
{

  // For now, we are not actually wrapping ComputePressure
  // To use EnzoComputePressure, we need to do some minor refactoring of it to
  // allow for optionally computing Pressure arrays specified in a Mapping.
  // This also requires making a modification to EnzoMethodGrackle's static
  // setup_grackle_fields method to also allow for specification of
  // relevant fields. Holding off on this for now

  if (grackle_variable_gamma_()){
    // we don't actually need to have grackle compute the pressure unless
    // gamma is allowed to vary
    ERROR("EnzoEOSIdeal::pressure_from_integrable",
	  "Not equipped to handle grackle and spatially variable gamma");
  }

  const bool idual = this->uses_dual_energy_formalism();
  const bool mag = (integrable_map.contains("bfield_x") ||
                    integrable_map.contains("bfield_y") ||
                    integrable_map.contains("bfield_z"));

  // rather than slicing out the unstaled regions, we may want use the full
  // array and adjust the iteration limits accordingly.

  EFlt3DArray density, vx, vy, vz, eint, etot, bx, by, bz;
  density = integrable_map.get("density", stale_depth);

  if (idual){
    eint = integrable_map.get("internal_energy", stale_depth);
  } else {
    etot = integrable_map.get("total_energy", stale_depth);
    vx = integrable_map.get("velocity_x", stale_depth);
    vy = integrable_map.get("velocity_y", stale_depth);
    vz = integrable_map.get("velocity_z", stale_depth);
    if (mag){
      bx = integrable_map.get("bfield_x", stale_depth);
      by = integrable_map.get("bfield_y", stale_depth);
      bz = integrable_map.get("bfield_z", stale_depth);
    }
  }

  CSlice unstaled(stale_depth,-stale_depth);
  EFlt3DArray p = pressure.subarray(unstaled, unstaled, unstaled);
  enzo_float gm1 = get_gamma() - 1.;

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	if (idual){
	  p(iz,iy,ix) = gm1 * density(iz,iy,ix) * eint(iz,iy,ix);
	} else {
          enzo_float v2 = (vx(iz,iy,ix) * vx(iz,iy,ix) +
			   vy(iz,iy,ix) * vy(iz,iy,ix) +
			   vz(iz,iy,ix) * vz(iz,iy,ix));
          enzo_float temp = (etot(iz,iy,ix) - 0.5 * v2) * density(iz,iy,ix);
          if (mag){
            enzo_float b2 = (bx(iz,iy,ix) * bx(iz,iy,ix) +
                             by(iz,iy,ix) * by(iz,iy,ix) +
                             bz(iz,iy,ix) * bz(iz,iy,ix));
            temp -= 0.5*b2;
          }
          p(iz,iy,ix) = gm1 * temp;
	}

      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoEOSIdeal::pressure_from_reconstructable
(EnzoEFltArrayMap &reconstructable, EFlt3DArray &pressure,
 int stale_depth) const
{
  // This is necessary since other equations of state may not include pressure
  // as a reconstructable quantity.

  // We will check if pressure_name the field containing pressure in
  // reconstructable_group are the same, if so then do nothing. Otherwise,
  // simply copy the values over.

  if (reconstructable.contains("pressure") &&
      (reconstructable.at("pressure").is_alias(pressure))){
    // The arrays are aliases of each other
    return;
  } else {
    EFlt3DArray old_p = reconstructable.get("pressure", stale_depth);

#ifdef DEBUG_MATCHING_ARRAY_SHAPES
    ASSERT6("EnzoEOSIdeal::pressure_from_reconstructable",
            ("The pressure array in reconstructable has shape (%d,%d,%d) "
             "while the array passed to this function has shape (%d,%d,%d). "
             "They should be the same."),
            old_p.shape(0), old_p.shape(1), old_p.shape(2),
            pressure.shape(0), pressure.shape(1), pressure.shape(2),
            ((old_p.shape(0) == pressure.shape(0)) &&
             (old_p.shape(1) == pressure.shape(1)) &&
             (old_p.shape(2) == pressure.shape(2))));
#endif

    for (int iz=0; iz< old_p.shape(0); iz++) {
      for (int iy=0; iy< old_p.shape(1); iy++) {
	for (int ix=0; ix< old_p.shape(2); ix++) {
	  pressure(iz,iy,ix) = old_p(iz,iy,ix);
	}
      }
    }
  }
}

//----------------------------------------------------------------------

// based on the enzo's hydro_rk implementation of synchronization (found in the
// Grid_UpdateMHD.C file)
void EnzoEOSIdeal::apply_floor_to_energy_and_sync
(EnzoEFltArrayMap &integrable_map, int stale_depth) const
{
  if (grackle_variable_gamma_()){
    ERROR("EnzoEOSIdeal::apply_floor_to_energy_and_sync",
	  "Not equipped to handle grackle and spatially variable gamma");
  }

  const bool idual = this->uses_dual_energy_formalism();
  const bool mag = (integrable_map.contains("bfield_x") ||
                    integrable_map.contains("bfield_y") ||
                    integrable_map.contains("bfield_z"));
  // in hydro_rk, eta was set equal to eta1 (it didn't use eta2 at all)
  const double eta = dual_energy_formalism_eta_;

  EFlt3DArray density, vx, vy, vz, etot, eint, bx, by, bz;
  density = integrable_map.get("density", stale_depth);
  vx = integrable_map.get("velocity_x", stale_depth);
  vy = integrable_map.get("velocity_y", stale_depth);
  vz = integrable_map.get("velocity_z", stale_depth);
  etot = integrable_map.get("total_energy", stale_depth);
  if (idual){
    eint = integrable_map.get("internal_energy", stale_depth);
  }
  if (mag){
    bx = integrable_map.get("bfield_x", stale_depth);
    by = integrable_map.get("bfield_y", stale_depth);
    bz = integrable_map.get("bfield_z", stale_depth);
  }

  float ggm1 = get_gamma()*(get_gamma() - 1.);
  enzo_float pressure_floor = get_pressure_floor();
  enzo_float inv_gm1 = 1./(get_gamma()-1.);

  // a requirement for an element of the internal energy field, cur_eint,
  // to be updated to the value computed from the total energy field, eint_1,
  // is that cur_eint > half_factor * cur_eint, where half_factor is 0.5. To
  // allow eta = 0, to specify that this update should always occur, we set
  // half_factor = 0 when eta = 0.
  const double half_factor = (eta != 0.) ? 0.5 : 0.;

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {

	enzo_float inv_rho = 1./density(iz,iy,ix);
	enzo_float eint_floor = pressure_floor*inv_gm1*inv_rho;

	enzo_float v2 = (vx(iz,iy,ix) * vx(iz,iy,ix) +
			 vy(iz,iy,ix) * vy(iz,iy,ix) +
			 vz(iz,iy,ix) * vz(iz,iy,ix));
        enzo_float non_thermal_e =  0.5*v2;
        enzo_float b2 = 0;
        if (mag){
          b2 = (bx(iz,iy,ix) * bx(iz,iy,ix) +
                by(iz,iy,ix) * by(iz,iy,ix) +
                bz(iz,iy,ix) * bz(iz,iy,ix));
          non_thermal_e += (0.5 * b2 *inv_rho);
        }

	if (idual){
	  enzo_float eint_1 = etot(iz,iy,ix) - non_thermal_e;
	  enzo_float cur_eint = eint(iz,iy,ix);

	  // compute cs^2 with estimate of eint from etot
	  // p = rho*(gamma-1)*eint
	  // cs^2 = gamma * p / rho = gamma*(gamma-1)*eint
	  enzo_float cs2_1 = std::fmax(0., ggm1*eint_1);

	  // half_factor = 0.5 when eta !=0. Otherwise it's 0.
	  if ( (cs2_1 > std::fmax(eta*v2, eta*b2*inv_rho)) &&
	       (eint_1 > half_factor*cur_eint) ){
	    cur_eint = eint_1;
	  }
	  cur_eint = EnzoEquationOfState::apply_floor(cur_eint, eint_floor);

	  eint(iz,iy,ix) = cur_eint;
	  etot(iz,iy,ix) = cur_eint + non_thermal_e;
	} else {

	  enzo_float etot_floor = eint_floor + non_thermal_e;
	  etot(iz,iy,ix) = EnzoEquationOfState::apply_floor(etot(iz,iy,ix),
							    etot_floor);
	}
      }
    }
  }
}
