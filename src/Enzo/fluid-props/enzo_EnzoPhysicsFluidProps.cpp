// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoPhysicsFluidProps.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-02-13
/// @brief    [\ref Enzo] Implementation of the EnzoPhysicsFluidProps class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

namespace { // define helper struct only used in this function
  struct NameVisitor {
    template<typename T> std::string operator()(T obj) { return T::name(); }
  };
}

enzo_float EnzoPhysicsFluidProps::gamma() const noexcept {

  // in c++17 this whole function simplifies down to:
  //
  // return eos_variant_.visit([](auto eos) -> enzo_float {
  //   using T = std::decay_t<decltype(eos)>;
  //   if constexpr (std::is_same<T,EnzoEOSIdeal>::value) { return eos.gamma; }
  //
  //   std::string name = T::name();
  //   ERROR1("EnzoPhysicsFluidProps::gamma()",
  //          "can't return gamma for EOS of type \"%s\". Either that type "
  //          "doesn't support it OR it's a new type & this function has not "
  //          "yet been updated to support it.", name.c_str());
  //   }
  // });

  if (eos_variant_.holds_alternative<EnzoEOSIdeal>()) {
    return eos_variant_.get<EnzoEOSIdeal>().gamma;
  }

  std::string name = eos_variant_.visit(NameVisitor());
  ERROR1("EnzoPhysicsFluidProps::gamma()",
         "can't return gamma for EOS of type \"%s\". Either that type "
         "doesn't support it OR it's a new type & this function has not "
         "yet been updated to support it.", name.c_str());
}

//----------------------------------------------------------------------

namespace { // define helper struct only used in this function
  struct IsBarotropicVisitor {
    template <typename T>
    bool operator()(T eos) const { return T::is_barotropic(); }
  };
}

bool EnzoPhysicsFluidProps::has_barotropic_eos() const noexcept {
  // in c++14 this function simplifies down to:
  // return eos_variant_.visit([](auto eos) { return eos.is_barotropic(); });
  return eos_variant_.visit(IsBarotropicVisitor());
}

//----------------------------------------------------------------------

void EnzoPhysicsFluidProps::primitive_from_integration
(const EnzoEFltArrayMap &integration_map, EnzoEFltArrayMap &primitive_map,
 const int stale_depth, const str_vec_t &passive_list,
 const bool ignore_grackle) const
{

  ASSERT("EnzoPhysicsFluidProps::primitive_from_integration",
         "the current implementation assumes a non-barotropic EOS",
         !this->has_barotropic_eos());

  // load shape ahead of time and declare them const for optimization purposes
  const CelloView<const enzo_float, 3> density = integration_map.at("density");
  const int mz = density.shape(0);
  const int my = density.shape(1);
  const int mx = density.shape(2);

  // This method was originally written as a method of the EOS class and since
  // then, we the EOS object didn't necessarily know what the integration
  // quantities are. For that reason, we take a somewhat exhaustive approach. 
  // - This operation could be more efficient if if were made a part of the
  //   Hydro/MHD Integrator
  // - alternatively, the EnzoPhysicsFluidProps should be able to infer this
  //   information in the future (if we pass this function hints about bfields
  //   and CRs)
  std::vector<std::string> quantity_list =
    EnzoCenteredFieldRegistry::get_registered_quantities(true, true);

  for (const std::string& key : quantity_list){
    if ((!integration_map.contains(key)) || (!primitive_map.contains(key))){
      continue;
    }

    const CelloView<const enzo_float, 3> integ_array = integration_map.at(key);
    const CelloView<enzo_float, 3> prim_array = primitive_map.at(key);

#ifdef DEBUG_MATCHING_ARRAY_SHAPES
    ASSERT6("EnzoPhysicsFluidProps::primitive_from_integration",
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
    const CelloView<const enzo_float, 3> cur_conserved
      = integration_map.at(key);
    const CelloView<enzo_float, 3> out_specific = primitive_map.at(key);

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

void EnzoPhysicsFluidProps::pressure_from_integration
(const EnzoEFltArrayMap &integration_map, const EFlt3DArray &pressure,
 const int stale_depth, const bool ignore_grackle) const
{

  const bool mhd = (integration_map.contains("bfield_x") ||
                    integration_map.contains("bfield_y") ||
                    integration_map.contains("bfield_z"));

  const bool idual = this->dual_energy_config().any_enabled();

  EnzoComputePressure::compute_pressure(EnzoFieldAdaptor(integration_map),
                                        pressure, mhd, idual, this->gamma(),
                                        stale_depth, ignore_grackle);
}

//----------------------------------------------------------------------

// based on the enzo's hydro_rk implementation of synchronization (found in the
// Grid_UpdateMHD.C file)
void EnzoPhysicsFluidProps::apply_floor_to_energy_and_sync
(EnzoEFltArrayMap &integration_map, const int stale_depth) const
{

  // This function's application of a floor isn't technically correct here for
  // a variable gamma.
  // - for (enzo::config()->method_grackle_chemistry->primordial_chemistry > 1)
  //   Grackle adjusts the "nominal gamma value" (usually ~ 5/3) based on the
  //   the relative abundance of molecular hydrogen & the specific internal
  //   energy (since the number of degrees of freedom depend on temperature)
  // - Grackle provides routines for calculating pressure and gamma given the
  //   mass_dens, eint, mass_dens_primordials, and mass_dens_H2.
  // - One could hypothetic invert the routine for pressure to acquire
  //   eint(mass_dens, pressure, mass_dens_primordials, and mass_dens_H2),
  //   but this is not presently available...
  //
  // The "correct" approach is to use the hypothetical eint function to compute
  // the local value of the internal energy floor for each cell using the
  // pressure_floor and the local values of mass_dens, mass_dens_primordials,
  // and mass_dens_H2.
  //
  // Since we don't have this hypothetical routine, we instead estimate the
  // local value of the internal energy floor for each cell using the pressure
  // floor, the "nominal gamma value", and the local mass_dens value, according
  // to eint = pressure / ((gamma - 1) * rho).
  // - This somewhat overestimates the true value of gamma.
  // - Thus, when you convert our eint_floor estimate back to pressure (with
  //   the Grackle routine), you'll recover a value smaller than the pressure
  //   floor.

  // this method does nothing for a barotropic eos
  if (this->has_barotropic_eos()) { return; }

  const EnzoDualEnergyConfig& de_config = this->dual_energy_config();

  const bool idual = de_config.any_enabled();

  if (idual && !de_config.modern_formulation()){
    ERROR("EnzoEOSIdeal::apply_floor_to_energy_and_sync",
          "The current implementation only works when the dual energy "
          "formalism is disabled or uses the \"modern formulation\"");
  }

  // retrieve the value of eta (if applicable)
  enzo_float tmp_eta = 0.0; // temporary variable
  de_config.modern_formulation(&tmp_eta);
  const double eta = tmp_eta;

  const bool mag = (integration_map.contains("bfield_x") ||
                    integration_map.contains("bfield_y") ||
                    integration_map.contains("bfield_z"));

  // historical context (in case we need to look back at the original code):
  // In enzo-dev's hydro_rk, eta was set equal to eta1 (it didn't ever use eta2)

  const EFlt3DArray etot = integration_map.at("total_energy");
  const EFlt3DArray eint =
    (idual) ? integration_map.at("internal_energy") : EFlt3DArray();

  using RdOnlyEFlt3DArray = CelloView<const enzo_float, 3>;
  const RdOnlyEFlt3DArray density = integration_map.at("density");
  const RdOnlyEFlt3DArray vx = integration_map.at("velocity_x");
  const RdOnlyEFlt3DArray vy = integration_map.at("velocity_y");
  const RdOnlyEFlt3DArray vz = integration_map.at("velocity_z");

  const RdOnlyEFlt3DArray bx = (mag) ?
    RdOnlyEFlt3DArray(integration_map.at("bfield_x")) : RdOnlyEFlt3DArray();
  const RdOnlyEFlt3DArray by = (mag) ?
    RdOnlyEFlt3DArray(integration_map.at("bfield_y")) : RdOnlyEFlt3DArray();
  const RdOnlyEFlt3DArray bz = (mag) ?
    RdOnlyEFlt3DArray(integration_map.at("bfield_z")) : RdOnlyEFlt3DArray();

  float ggm1 = this->gamma()*(this->gamma() - 1.);
  enzo_float pressure_floor = this->fluid_floor_config().pressure();
  enzo_float inv_gm1 = 1./(this->gamma()-1.);

  // a requirement for an element of the internal energy field, cur_eint,
  // to be updated to the value computed from the total energy field, eint_1,
  // is that cur_eint > half_factor * cur_eint, where half_factor is 0.5. To
  // allow eta = 0, to specify that this update should always occur, we set
  // half_factor = 0 when eta = 0.
  const double half_factor = (eta != 0.) ? 0.5 : 0.;

  for (int iz = stale_depth; iz < (density.shape(0) - stale_depth); iz++) {
    for (int iy = stale_depth; iy < (density.shape(1) - stale_depth); iy++) {
      for (int ix = stale_depth; ix < (density.shape(2) - stale_depth); ix++) {

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
	  cur_eint = enzo_utils::apply_floor(cur_eint, eint_floor);

	  eint(iz,iy,ix) = cur_eint;
	  etot(iz,iy,ix) = cur_eint + non_thermal_e;
	} else {

	  enzo_float etot_floor = eint_floor + non_thermal_e;
	  etot(iz,iy,ix) = enzo_utils::apply_floor(etot(iz,iy,ix), etot_floor);
	}
      }
    }
  }
}
