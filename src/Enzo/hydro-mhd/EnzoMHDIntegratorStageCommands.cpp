// See LICENSE_CELLO file for license and copyright information

/// @file     EnzoMHDIntegratorStageCommands.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Sun May 28 2023
/// @brief    [\ref Enzo] Implementation of EnzoMHDIntegratorStageCommands

#include "Enzo/hydro-mhd/EnzoMHDIntegratorStageCommands.hpp"
#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/hydro-mhd/hydro-mhd.hpp"

#include <algorithm> // std::transform
#include <cctype> // tolower

//----------------------------------------------------------------------

EnzoMHDIntegratorStageCommands::EnzoMHDIntegratorStageCommands
(const EnzoMHDIntegratorStageArgPack& args)
  : riemann_solver_(nullptr),
    reconstructors_(),
    integration_quan_updater_(nullptr),
    mhd_choice_(EnzoMHDIntegratorStageCommands::parse_bfield_choice_
                (args.mhd_choice))
{
  // check compatability with EnzoPhysicsFluidProps
  EnzoPhysicsFluidProps* fluid_props = enzo::fluid_props();
  ASSERT("EnzoMethodMHDVlct::EnzoMethodMHDVlct",
         "can't currently handle the case with a non-ideal EOS",
         fluid_props->eos_variant().holds_alternative<EnzoEOSIdeal>());
  const EnzoDualEnergyConfig& de_config = fluid_props->dual_energy_config();
  ASSERT("EnzoMethodMHDVlct::EnzoMethodMHDVlct",
         "selected formulation of dual energy formalism is incompatible",
         de_config.is_disabled() || de_config.modern_formulation());
  const EnzoFluidFloorConfig& fluid_floor_config
    = fluid_props->fluid_floor_config();
  ASSERT("EnzoMethodMHDVlct::EnzoMethodMHDVlct",
         "density and pressure floors must be defined",
         fluid_floor_config.has_density_floor() &
         fluid_floor_config.has_pressure_floor());

  // Initialize the Riemann Solver
  riemann_solver_ = EnzoRiemann::construct_riemann
    ({args.rsolver, mhd_choice_ != bfield_choice::no_bfield,
      de_config.any_enabled()});

  // determine integration and primitive field list
  str_vec_t integration_field_list
    = riemann_solver_->integration_quantity_keys();
  str_vec_t primitive_field_list
    = riemann_solver_->primitive_quantity_keys();

  // Initialize the remaining component objects
  const enzo_float theta_limiter = static_cast<enzo_float>(args.theta_limiter);
  for (const std::string& name : args.recon_names) {
    std::unique_ptr<EnzoReconstructor> temp
      (EnzoReconstructor::construct_reconstructor(primitive_field_list, name,
                                                  theta_limiter));
    reconstructors_.push_back(std::move(temp));
  }

  integration_quan_updater_ =
    new EnzoIntegrationQuanUpdate(integration_field_list, true);
}

//----------------------------------------------------------------------

EnzoMHDIntegratorStageCommands::~EnzoMHDIntegratorStageCommands()
{
  delete riemann_solver_;
  delete integration_quan_updater_;
}

//----------------------------------------------------------------------

EnzoMHDIntegratorStageCommands::bfield_choice
EnzoMHDIntegratorStageCommands::parse_bfield_choice_
(std::string choice) noexcept
{
  std::string formatted(choice.size(), ' ');
  std::transform(choice.begin(), choice.end(), formatted.begin(),
		 ::tolower);
  if (formatted == std::string("no_bfield")){
    return bfield_choice::no_bfield;
  } else if (formatted == std::string("unsafe_constant_uniform")){
    ERROR("EnzoMHDIntegratorStageCommands::parse_bfield_choice_",
          "constant_uniform is primarilly for debugging purposes. DON'T use "
          "for science runs (things can break).");
    return bfield_choice::unsafe_const_uniform;
  } else if (formatted == std::string("constrained_transport")){
    return bfield_choice::constrained_transport;
  } else {
    ERROR("EnzoMHDIntegratorStageCommands::parse_bfield_choice_",
          "Unrecognized choice. Known options include \"no_bfield\" and "
          "\"constrained_transport\"");
    return bfield_choice::no_bfield;
  }
}

//----------------------------------------------------------------------

void EnzoMHDIntegratorStageCommands::compute_update_stage
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
 int stale_depth,
 const std::array<enzo_float,3> cell_widths_xyz
 ) const noexcept
{
  check_valid_stage_index_(static_cast<int>(stage_index));

  EnzoPhysicsFluidProps* fluid_props = enzo::fluid_props();

  EnzoReconstructor *reconstructor = reconstructors_[stage_index].get();

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
  }

  // Update the integration quantities (includes flux divergence and source
  // terms). This currently needs to happen after updating the
  // cell-centered B-field so that the pressure floor can be applied to the
  // total energy (and if necessary the total energy can be synchronized
  // with the internal energy)
  integration_quan_updater_->update_quantities
    (tstep_begin_integration_map, dUcons_map, out_integration_map,
     stale_depth, passive_list);
}

//----------------------------------------------------------------------

void EnzoMHDIntegratorStageCommands::compute_flux_
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

void EnzoMHDIntegratorStageCommands::compute_source_terms_
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

double EnzoMHDIntegratorStageCommands::timestep
(EnzoEFltArrayMap integration_map,
 CelloView<const enzo_float, 3> pressure,
 const double dx, const double dy, const double dz) const noexcept
{

  using RdOnlyEFltView = CelloView<const enzo_float, 3>;
  const RdOnlyEFltView density = integration_map.at("density");

  const RdOnlyEFltView velocity_x = integration_map.at("velocity_x");
  const RdOnlyEFltView velocity_y = integration_map.at("velocity_y");
  const RdOnlyEFltView velocity_z = integration_map.at("velocity_z");

  // this will raise an error if not Ideal EOS
  const EnzoEOSIdeal eos
    = enzo::fluid_props()->eos_variant().get<EnzoEOSIdeal>();

  // timestep is the minimum of dr_i/(abs(v_i)+signal_speed) for all dimensions
  //   - dr_i and v_i are the the width of the cell and velocity along
  //     dimension i.
  //   - signal_speed is sound speed without Bfields and fast-magnetosonic
  //     speed with Bfields

  // initialize returned value
  double dtBaryons = std::numeric_limits<double>::max();

  if (this->is_pure_hydro()) {
    auto loop_body = [=, &dtBaryons](int iz, int iy, int ix)
      {
        double cs = (double) eos.sound_speed(density(iz,iy,ix),
                                             pressure(iz,iy,ix));
        double local_dt = enzo_utils::min<double>
          (dx/(std::fabs((double) velocity_x(iz,iy,ix)) + cs),
           dy/(std::fabs((double) velocity_y(iz,iy,ix)) + cs),
           dz/(std::fabs((double) velocity_z(iz,iy,ix)) + cs));
        dtBaryons = std::min(dtBaryons, local_dt);
      };
    enzo_utils::exec_loop(density.shape(0), density.shape(1), density.shape(2),
                          0, loop_body);

  } else {
    const RdOnlyEFltView bfieldc_x = integration_map.at("bfield_x");
    const RdOnlyEFltView bfieldc_y = integration_map.at("bfield_y");
    const RdOnlyEFltView bfieldc_z = integration_map.at("bfield_z");

    auto loop_body = [=, &dtBaryons](int iz, int iy, int ix)
      {
        // if the bfield is 0 at any given point, the fast magnetosonic speed
        // correctly reduces to the sound speed.
        //
        // We follow the convention of using the maximum value of the fast
        // magnetosonic speed:     cfast = sqrt(va^2+cs^2)
        double cfast = (double) eos.fast_magnetosonic_speed<0>
          (density(iz,iy,ix), pressure(iz,iy,ix),
           bfieldc_x(iz,iy,ix), bfieldc_y(iz,iy,ix), bfieldc_z(iz,iy,ix));

        double local_dt = enzo_utils::min<double>
           (dx/(std::fabs((double) velocity_x(iz,iy,ix)) + cfast),
            dy/(std::fabs((double) velocity_y(iz,iy,ix)) + cfast),
            dz/(std::fabs((double) velocity_z(iz,iy,ix)) + cfast));
        dtBaryons = std::min(dtBaryons, local_dt);
      };
    enzo_utils::exec_loop(density.shape(0), density.shape(1), density.shape(2),
                          0, loop_body);
  }

  return dtBaryons; // courant factor is handled separately!
}
