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

void EnzoMHDVlctIntegrator::compute_update_stage
(EnzoEFltArrayMap tstep_begin_integration_map,
 EnzoEFltArrayMap cur_stage_integration_map,
 EnzoEFltArrayMap out_integration_map,
 EnzoEFltArrayMap primitive_map,
 EnzoEFltArrayMap priml_map, EnzoEFltArrayMap primr_map,
 std::array<EnzoEFltArrayMap, 3> flux_maps_xyz,
 EnzoEFltArrayMap dUcons_map,
 const EnzoEFltArrayMap accel_map,
 EFlt3DArray interface_vel_arr,
 const str_vec_t& passive_list,
 EnzoBfieldMethod* bfield_method_,
 unsigned short stage_index,
 double cur_dt,
 int& stale_depth, // to do: come up with a way to more directly communicate
                   // changes in stale_depth
 const std::array<enzo_float,3> cell_widths_xyz
 )
{
  ASSERT("EnzoMHDVlctIntegrator::compute_update_stage",
         "stage_index must satisfy 0 <= stage_index < 2", stage_index < 2);

  EnzoPhysicsFluidProps* fluid_props = enzo::fluid_props();

  EnzoReconstructor *reconstructor =
    (stage_index == 0) ? half_dt_recon_ : full_dt_recon_;

  // set all elements of the arrays in dUcons_map to 0 (throughout the rest
  // of the current loop, flux divergence and source terms will be
  // accumulated in these arrays)
  integration_quan_updater_->clear_dUcons_map(dUcons_map, 0., passive_list);

  // Compute the primitive quantities from the integration quantites
  // This basically copies all quantities that are both an integration
  // quantity and a primitive and converts the passsive scalars from
  // conserved-form to specific-form (i.e. from density to mass fraction).
  // For a non-barotropic gas, this also computes pressure
  // - for consistency with the Ppm solver, we explicitly avoid the Grackle
  //   routine. This is only meaningful when grackle models molecular
  //   hydrogen (which modifes the adiabtic index)
  const bool ignore_grackle = true;
  fluid_props->primitive_from_integration(cur_stage_integration_map,
                                          primitive_map,
                                          stale_depth, passive_list,
                                          ignore_grackle);

  // Compute flux along each dimension
  for (int dim = 0; dim < 3; dim++){
    // trim the shape of priml_map and primr_map (they're bigger than
    // necessary so that they can be reused for each dim).
    CSlice x_slc = (dim == 0) ? CSlice(0,-1) : CSlice(0, nullptr);
    CSlice y_slc = (dim == 1) ? CSlice(0,-1) : CSlice(0, nullptr);
    CSlice z_slc = (dim == 2) ? CSlice(0,-1) : CSlice(0, nullptr);

    EnzoEFltArrayMap pl_map = priml_map.subarray_map(z_slc, y_slc, x_slc);
    EnzoEFltArrayMap pr_map = primr_map.subarray_map(z_slc, y_slc, x_slc);

    EFlt3DArray *interface_vel_arr_ptr, sliced_interface_vel_arr;
    if (fluid_props->dual_energy_config().any_enabled()){
      // when using dual energy formalism, trim the trim scratch-array for
      // storing interface velocity values (computed by the Riemann Solver).
      // This is used in the calculation of the internal energy source
      // term). As with priml_map and primr_map, the array is bigger than
      // necessary so it can be reused for each axis
      sliced_interface_vel_arr =
        interface_vel_arr.subarray(z_slc, y_slc, x_slc);
      interface_vel_arr_ptr = &sliced_interface_vel_arr;
    } else {
      // no scratch-space was allocated, so we just pass a nullptr
      interface_vel_arr_ptr = nullptr;
    }

    compute_flux_(dim, cur_dt, cell_widths_xyz[dim], primitive_map,
                  pl_map, pr_map, flux_maps_xyz[dim], dUcons_map,
                  interface_vel_arr_ptr, *reconstructor, bfield_method_,
                  stale_depth, passive_list);
  }

  // increment the stale_depth
  stale_depth+=reconstructor->immediate_staling_rate();

  // Compute the source terms (use them to update dUcons_group)
  compute_source_terms_(cur_dt, (stage_index == 1), tstep_begin_integration_map,
                        primitive_map, accel_map, dUcons_map, stale_depth);

  // Update Bfields
  if (bfield_method_ != nullptr) {
    bfield_method_->update_all_bfield_components(cur_stage_integration_map,
                                                 flux_maps_xyz[0],
                                                 flux_maps_xyz[1],
                                                 flux_maps_xyz[2],
                                                 out_integration_map,
                                                 cur_dt,
                                                 stale_depth);
    bfield_method_->increment_partial_timestep();
  }

  // Update the integration quantities (includes flux divergence and source
  // terms). This currently needs to happen after updating the
  // cell-centered B-field so that the pressure floor can be applied to the
  // total energy (and if necessary the total energy can be synchronized
  // with the internal energy)
  integration_quan_updater_->update_quantities
    (tstep_begin_integration_map, dUcons_map, out_integration_map,
     stale_depth, passive_list);

  // increment stale_depth since the inner values have been updated
  // but the outer values have not
  stale_depth+=reconstructor->delayed_staling_rate();
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

void EnzoMHDVlctIntegrator::compute_source_terms_
(const double cur_dt, const bool full_timestep,
 const EnzoEFltArrayMap &orig_integration_map,
 const EnzoEFltArrayMap &primitive_map,
 const EnzoEFltArrayMap &accel_map,
 EnzoEFltArrayMap &dUcons_map,  const int stale_depth) const noexcept
{
  // As we add more source terms, we may want to store them in a std::vector
  // instead of manually invoking them

  // add any source-terms that are used for partial and full timesteps
  // (there aren't any right now...)

  // add any source-terms that are only included for full-timestep
  if (full_timestep & (accel_map.size() != 0)){
    // include gravity source terms.
    //
    // The inclusion of these terms are not based on any external paper. Thus,
    // we include the following bullets to explain our thought process:
    // - to be safe, we will use the density & velocity values from the start
    //   of the timestep to compute the source terms.
    // - we should reconsider this choice if we later decide to include the
    //   source term for both the partial & full timesteps.
    // - to include the source term during the partial & full timesteps, we
    //   would probably want to recompute the gravitational potential and
    //   acceleration fields at the partial timestep
    EnzoSourceGravity gravity_source;
    gravity_source.calculate_source(cur_dt, orig_integration_map, dUcons_map,
                                    accel_map, stale_depth);
  }
}

//----------------------------------------------------------------------
