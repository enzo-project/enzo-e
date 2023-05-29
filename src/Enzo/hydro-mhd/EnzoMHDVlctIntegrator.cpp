// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMHDVlctIntegrator.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 28 2023
/// @brief    [\ref Enzo] Implementation of EnzoMHDVlctIntegrator

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"

// place this after #include "enzo.hpp"
#include "EnzoMHDVlctIntegrator.hpp"

//----------------------------------------------------------------------

EnzoMHDVlctIntegrator::~EnzoMHDVlctIntegrator()
{
  delete half_dt_recon_;
  delete full_dt_recon_;
  delete riemann_solver_;
  delete integration_quan_updater_;
}

//----------------------------------------------------------------------

void EnzoMHDVlctIntegrator::compute_flux_
(const int dim, const double cur_dt, const enzo_float cell_width,
 EnzoEFltArrayMap &primitive_map,
 EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
 EnzoEFltArrayMap &flux_map, EnzoEFltArrayMap &dUcons_map,
 const EFlt3DArray* const interface_velocity_arr_ptr,
 EnzoReconstructor &reconstructor, EnzoBfieldMethod *bfield_method,
 const int stale_depth, const str_vec_t& passive_list) const noexcept
{

  // First, reconstruct the left and right interface values
  reconstructor.reconstruct_interface(primitive_map, priml_map, primr_map,
				      dim, stale_depth, passive_list);

  // We temporarily increment the stale_depth for the rest of this calculation
  // here. We can't fully increment otherwise it will screw up the
  // reconstruction along other dimensions
  int cur_stale_depth = stale_depth + reconstructor.immediate_staling_rate();

  // Need to set the component of reconstructed B-field along dim, equal to
  // the corresponding longitudinal component of the B-field tracked at cell
  // interfaces
  if (bfield_method != nullptr) {
    bfield_method->correct_reconstructed_bfield(priml_map, primr_map,
                                                dim, cur_stale_depth);
  }

  // Next, compute the fluxes
  riemann_solver_->solve(priml_map, primr_map, flux_map, dim,
                         cur_stale_depth, passive_list,
                         interface_velocity_arr_ptr);

  // Compute the changes in the conserved form of the integration quantities
  // from the fluxes and use these values to update dUcons_map (which is used
  // to accumulate the total change in these quantities over the current
  // [partial] timestep)
  integration_quan_updater_->accumulate_flux_component(dim, cur_dt, cell_width,
                                                       flux_map, dUcons_map,
                                                       cur_stale_depth,
                                                       passive_list);

  // if using dual energy formalism, compute the component of the internal
  // energy source term for this dim (and update dUcons_map).
  if (enzo::fluid_props()->dual_energy_config().any_enabled()){
    EnzoSourceInternalEnergy eint_src;
    eint_src.calculate_source(dim, cur_dt, cell_width, primitive_map,
                              dUcons_map, *interface_velocity_arr_ptr,
                              cur_stale_depth);
  }

  // Finally, have bfield_method record the upwind direction (for handling CT)
  if (bfield_method != nullptr){
    bfield_method->identify_upwind(flux_map, dim, cur_stale_depth);
  }
}

//----------------------------------------------------------------------
