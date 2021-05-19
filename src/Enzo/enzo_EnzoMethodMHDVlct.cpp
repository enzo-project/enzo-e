// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodMHDVlct.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri June 14 2019
/// @brief    [\ref Enzo] Implementation of the EnzoMethodMHDVlct class

#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"
#include <algorithm>    // std::copy

//----------------------------------------------------------------------

static void build_field_l_(const std::vector<std::string> &quantity_l,
                           std::vector<std::string> &field_l)
{
  for (const std::string& quantity : quantity_l){
    bool success, is_vector;
    success = EnzoCenteredFieldRegistry::quantity_properties (quantity,
                                                              &is_vector);

    ASSERT1("build_field_l_",
            ("\"%s\" is not registered in EnzoCenteredFieldRegistry"),
            quantity.c_str(), success);

    if (is_vector){
      field_l.push_back(quantity + "_x");
      field_l.push_back(quantity + "_y");
      field_l.push_back(quantity + "_z");
    } else {
      field_l.push_back(quantity);
    }
  }

  FieldDescr * field_descr = cello::field_descr();
  for (const std::string& field : field_l){
    ASSERT1("EnzoMethodMHDVlct", "\"%s\" must be a permanent field",
	    field.c_str(), field_descr->is_field(field));
  }
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::EnzoMethodMHDVlct (std::string rsolver,
				      std::string half_recon_name,
				      std::string full_recon_name,
				      double gamma, double theta_limiter,
				      double density_floor,
				      double pressure_floor,
				      std::string mhd_choice,
				      bool dual_energy_formalism,
				      double dual_energy_formalism_eta)
  : Method()
{
  // Initialize equation of state (check the validity of quantity floors)
  EnzoEquationOfState::check_floor(density_floor);
  EnzoEquationOfState::check_floor(pressure_floor);
  eos_ = new EnzoEOSIdeal(gamma, density_floor, pressure_floor,
			  dual_energy_formalism, dual_energy_formalism_eta);

  // Determine whether magnetic fields are to be used
  mhd_choice_ = parse_bfield_choice_(mhd_choice);

  // determine integrable and reconstructable quantities
  std::vector<std::string> integrable_quantities, reconstructable_quantities;
  EnzoMethodMHDVlct::determine_quantities_(eos_, mhd_choice_,
                                           integrable_quantities,
                                           reconstructable_quantities);

  // Initialize the remaining component objects
  half_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (reconstructable_quantities, half_recon_name, (enzo_float)theta_limiter);
  full_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (reconstructable_quantities, full_recon_name, (enzo_float)theta_limiter);
  riemann_solver_ = EnzoRiemann::construct_riemann(integrable_quantities,
                                                   rsolver);
  integrable_updater_ = new EnzoIntegrableUpdate(integrable_quantities,
                                                 true);

  // Determine the list of integrable fields and reconstructable fields and
  // then make sure all required fields are defined
  build_field_l_(integrable_quantities, integrable_field_list_);
  build_field_l_(reconstructable_quantities, reconstructable_field_list_);

  // make sure "pressure" is defined (it's needed to compute the timestep)
  FieldDescr * field_descr = cello::field_descr();
  ASSERT("EnzoMethodMHDVlct", "\"pressure\" must be a permanent field",
	 field_descr->is_field("pressure"));

  if (mhd_choice_ == bfield_choice::constrained_transport) {
    bfield_method_ = new EnzoBfieldMethodCT(2);
    bfield_method_->check_required_fields();
  } else {
    bfield_method_ = nullptr;
  }

  // Finally, initialize the default Refresh object
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  // Need to refresh all fields because the fields holding passively advected
  // scalars won't necessarily be known until after all Methods have been
  // constructed and all intializers have been executed
  refresh->add_all_fields();
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::bfield_choice EnzoMethodMHDVlct::parse_bfield_choice_
(std::string choice) const noexcept
{
  std::string formatted(choice.size(), ' ');
  std::transform(choice.begin(), choice.end(), formatted.begin(),
		 ::tolower);
  if (formatted == std::string("no_bfield")){
    return bfield_choice::no_bfield;
  } else if (formatted == std::string("unsafe_constant_uniform")){
    ERROR("EnzoMethodMHDVlct::parse_bfield_choice_",
          "constant_uniform is primarilly for debugging purposes. DON'T use "
          "for science runs (things can break).");
    return bfield_choice::unsafe_const_uniform;
  } else if (formatted == std::string("constrained_transport")){
    return bfield_choice::constrained_transport;
  } else {
    ERROR("EnzoMethodMHDVlct::parse_bfield_choice_",
          "Unrecognized choice. Known options include \"no_bfield\" and "
          "\"constrained_transport\"");
    return bfield_choice::no_bfield;
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::determine_quantities_
(const EnzoEquationOfState *eos, EnzoMethodMHDVlct::bfield_choice mhd_choice,
 std::vector<std::string> &integrable_quantities,
 std::vector<std::string> &reconstructable_quantities) noexcept
{
#ifdef CONFIG_USE_GRACKLE
  if (enzo::config()->method_grackle_use_grackle){
    // we can remove the following once EnzoMethodGrackle no longer requires
    // the internal_energy to be a permanent field
    ASSERT("EnzoMethodMHDVlct::determine_quantities_",
           ("Grackle cannot currently be used alongside this integrator "
            "unless the dual-energy formalism is in use"),
           eos->uses_dual_energy_formalism());
  }
#endif /* CONFIG_USE_GRACKLE */

  std::string common[] {"density", "velocity"};
  for (std::string quantity : common){
    integrable_quantities.push_back(quantity);
    reconstructable_quantities.push_back(quantity);
  }

  if (mhd_choice != bfield_choice::no_bfield){
    integrable_quantities.push_back("bfield");
    reconstructable_quantities.push_back("bfield");
  }

  if (eos->is_barotropic()){
    ERROR("EnzoMethodMHDVlct::determine_quantities_",
	  "Not presently equipped to handle barotropic equations of state.");
  } else {
    // add specific total energy to integrable quantities
    integrable_quantities.push_back("total_energy");
    // add pressure to reconstructable quantities
    reconstructable_quantities.push_back("pressure");
    // add specific internal energy to integrable quantities (if using the dual
    // engery formalism)
    if (eos->uses_dual_energy_formalism()){
      integrable_quantities.push_back("internal_energy");
    }
  }
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::~EnzoMethodMHDVlct()
{
  delete eos_;
  delete half_dt_recon_;
  delete full_dt_recon_;
  delete riemann_solver_;
  delete integrable_updater_;
  if (bfield_method_ != nullptr){
    delete bfield_method_;
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p|eos_;
  p|half_dt_recon_;
  p|full_dt_recon_;
  p|riemann_solver_;
  p|integrable_updater_;
  p|mhd_choice_;
  p|bfield_method_;
  p|integrable_field_list_;
  p|reconstructable_field_list_;
  p|lazy_passive_list_;
}

//----------------------------------------------------------------------

// Returns the unique members of a combination of 2 vectors
std::vector<std::string> unique_combination_(const std::vector<std::string> &a,
					     const std::vector<std::string> &b)
{
  std::vector<std::string> out;
  // copy elements from a
  for (std::size_t i = 0; i<a.size(); i++){
    std::string name = a[i];
    if (std::find(out.begin(), out.end(), name) == out.end()){
      out.push_back(name);
    }
  }
  // copy elements from b
  for (std::size_t i = 0; i<b.size(); i++){
    std::string name = b[i];
    if (std::find(out.begin(), out.end(), name) == out.end()){
      out.push_back(name);
    }
  }
  return out;
}

//----------------------------------------------------------------------

EnzoEFltArrayMap EnzoMethodMHDVlct::nonpassive_primitive_map_(Block * block)
  const noexcept
{
  EnzoEFltArrayMap map("primitive");
  std::vector<std::string> all_prim_fields =
    unique_combination_(integrable_field_list_, reconstructable_field_list_);
  EnzoFieldArrayFactory array_factory(block,0);
  for (const std::string& field_name : all_prim_fields){
    ASSERT1("EnzoMethodMHDVlct::nonpassive_primitive_map_",
            "EnzoEFltArrayMap can't hold more than one key called \"%s\"",
            field_name.c_str(), !map.contains(field_name));
    map[field_name] = array_factory.from_name(field_name);
  }
  return map;
}

//----------------------------------------------------------------------

EnzoEFltArrayMap EnzoMethodMHDVlct::conserved_passive_scalar_map_
(Block * block) const noexcept
{
  EnzoEFltArrayMap map("conserved_passive_scalar");
  std::shared_ptr<const str_vec_t> passive_list= lazy_passive_list_.get_list();
  EnzoFieldArrayFactory array_factory(block,0);
  for (const std::string& field_name : (*passive_list)){
    ASSERT1("EnzoMethodMHDVlct::conserved_passive_scalar_map_",
            "EnzoEFltArrayMap can't hold more than one key called \"%s\"",
            field_name.c_str(), !map.contains(field_name));
    map[field_name] = array_factory.from_name(field_name);
  }
  return map;
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute ( Block * block) throw()
{
  if (block->is_leaf()) {
    // Check that the mesh size and ghost depths are appropriate
    check_mesh_and_ghost_size_(block);

    // declaring Maps of arrays and stand-alone arrays that wrap existing
    // fields and/or serve as scratch space.

    // map that holds arrays wrapping the Cello Fields holding each of the
    // primitive quantities. Additionally, this also includes temporary arrays
    // used to hold the specific form of the passive scalar
    EnzoEFltArrayMap primitive_map; // this will be overwritten

    // map used for storing primitive values at the half time-step. This
    // includes key,array pairs for each entry in primitive_map (there should
    // be no aliased fields shared between maps)
    EnzoEFltArrayMap temp_primitive_map("temp_primitive");

    // map holding the arrays wrapping each fields corresponding to a passively
    // advected scalar. The scalar is in conserved form.
    EnzoEFltArrayMap conserved_passive_scalar_map; // this will be overwritten

    // holds left and right reconstructed primitives (scratch-space)
    EnzoEFltArrayMap priml_map("priml");
    EnzoEFltArrayMap primr_map("primr");

    // Arrays used to store the pressure computed from the reconstructed left
    // and right primitives
    // Note: in the case of adiabatic fluids, pressure is a reconstructable
    //       quantity and entries are included for it in priml_map and
    //       primr_map. In that case, pressure_l and pressure_r are aliases of
    //       those arrays.
    EFlt3DArray pressure_l, pressure_r;

    // maps used to store fluxes (in the future, these will wrap FluxData
    // entries)
    EnzoEFltArrayMap xflux_map("xflux");
    EnzoEFltArrayMap yflux_map("yflux");
    EnzoEFltArrayMap zflux_map("zflux");

    // map of arrays  used to accumulate the changes to the conserved forms of
    // the integrable quantities and passively advected scalars. In other
    // words, at the start of the (partial) timestep, the fields are all set to
    // zero and are used to accumulate the flux divergence and source terms. If
    // CT is used, it won't have space to store changes in the magnetic fields.
    EnzoEFltArrayMap dUcons_map("dUcons");

    setup_arrays_(block, primitive_map, temp_primitive_map,
                  conserved_passive_scalar_map, priml_map, primr_map,
		  pressure_l, pressure_r, xflux_map, yflux_map, zflux_map,
                  dUcons_map);

    // Setup a pointer to an array that used to store interface velocity fields
    // from computed by the Riemann Solver (to use in the calculation of the
    // internal energy source term). If the dual energy formalism is not in
    // use, don't actually allocate the array and set the pointer to NULL.
    EFlt3DArray interface_velocity_arr, *interface_velocity_arr_ptr;
    if (eos_->uses_dual_energy_formalism()){
      EFlt3DArray density = primitive_map.at("density");
      interface_velocity_arr = EFlt3DArray(density.shape(0), density.shape(1),
                                           density.shape(2));
      interface_velocity_arr_ptr = &interface_velocity_arr;
    } else {
      interface_velocity_arr_ptr = nullptr;
    }

    // allocate constrained transport object
    if (bfield_method_ != nullptr) {
      bfield_method_->register_target_block(block);
    }

    const enzo_float* const cell_widths = enzo::block(block)->CellWidth;

    double dt = block->dt();

    // stale_depth indicates the number of field entries from the outermost
    // field value that the region including "stale" values (need to be
    // refreshed) extends over.
    int stale_depth = 0;

    // convert the passive scalars from conserved form to specific form
    // (outside the integrator, they are treated like conserved densities)
    compute_specific_passive_scalars_(*(lazy_passive_list_.get_list()),
                                      primitive_map["density"],
                                      conserved_passive_scalar_map,
                                      primitive_map, stale_depth);

    // repeat the following loop twice (for half time-step and full time-step)

    for (int i=0;i<2;i++){
      double cur_dt = (i == 0) ? dt/2. : dt;
      EnzoEFltArrayMap& cur_integrable_map =
        (i == 0) ? primitive_map      : temp_primitive_map;
      EnzoEFltArrayMap& out_integrable_map =
        (i == 0) ? temp_primitive_map :      primitive_map;
      // For the purposes of making the calculation procedure slightly more
      // more explicit, we distinguish between cur_integrable_group and
      // cur_reconstructable_group. Due to the high level of overlap between
      // these, they are simply aliases of the same underlying Grouping that
      // holds groups for both of them
      EnzoEFltArrayMap& cur_reconstructable_map = cur_integrable_map;

      EnzoReconstructor *reconstructor;

      if (i == 0){
        reconstructor = half_dt_recon_;
      } else {
        reconstructor = full_dt_recon_;

        // After the fluxes were added to the passive scalar in the first half
        // timestep, the values were stored in conserved form in the fields
        // held by conserved_passive_scalar_map.
        // Need to convert them to specific form
        compute_specific_passive_scalars_(*(lazy_passive_list_.get_list()),
                                          cur_integrable_map["density"],
                                          conserved_passive_scalar_map,
                                          cur_integrable_map, stale_depth);
      }

      // set all elements of the arrays in dUcons_map to 0 (throughout the rest
      // of the current loop, flux divergence and source terms will be
      // accumulated in these arrays)
      integrable_updater_->clear_dUcons_map(dUcons_map, 0.,
                                            *(lazy_passive_list_.get_list()));

      // Compute the reconstructable quantities from the integrable quantites
      // Although cur_integrable_map holds the passive scalars in integrable
      // form, the conserved form of the values is required in case Grackle is
      // being used.
      //
      // Note: cur_integrable_map and cur_reconstructable_map are aliases of
      // the same map since there is such a large degree of overlap between
      // reconstructable and integrable quantities
      //
      // For a barotropic gas, the following nominally does nothing
      // For a non-barotropic gas, the following nominally computes pressure
      eos_->reconstructable_from_integrable(cur_integrable_map,
                                            cur_reconstructable_map,
                                            conserved_passive_scalar_map,
                                            stale_depth,
                                            *(lazy_passive_list_.get_list()));

      // Compute flux along each dimension
      compute_flux_(0, cur_dt, cell_widths[0], cur_reconstructable_map,
                    priml_map, primr_map, pressure_l, pressure_r,
                    xflux_map, dUcons_map, interface_velocity_arr_ptr,
                    *reconstructor, bfield_method_, stale_depth,
                    *(lazy_passive_list_.get_list()));
      compute_flux_(1, cur_dt, cell_widths[1], cur_reconstructable_map,
                    priml_map, primr_map, pressure_l, pressure_r,
                    yflux_map, dUcons_map, interface_velocity_arr_ptr,
                    *reconstructor, bfield_method_, stale_depth,
                    *(lazy_passive_list_.get_list()));
      compute_flux_(2, cur_dt, cell_widths[2], cur_reconstructable_map,
                    priml_map, primr_map, pressure_l, pressure_r,
                    zflux_map, dUcons_map, interface_velocity_arr_ptr,
                    *reconstructor, bfield_method_, stale_depth,
                    *(lazy_passive_list_.get_list()));

      // increment the stale_depth
      stale_depth+=reconstructor->immediate_staling_rate();

      // This is where source terms should be computed (added to dUcons_group)

      // Update Bfields
      if (bfield_method_ != nullptr) {
        bfield_method_->update_all_bfield_components(cur_integrable_map,
                                                     xflux_map, yflux_map,
                                                     zflux_map,
                                                     out_integrable_map, cur_dt,
                                                     stale_depth);
        bfield_method_->increment_partial_timestep();
      }

      // Update quantities (includes flux divergence and source terms) 
      // This needs to happen after updating the cell-centered B-field so that
      // the pressure floor can be applied to the total energy (and if
      // necessary the total energy can be synchronized with internal energy)
      //
      // Note: updated passive scalars are NOT saved in out_integrable_group in
      //     specific form. Instead they are saved in conserved_passive_scalars
      //     in conserved form.
      integrable_updater_->update_quantities
        (primitive_map, dUcons_map, out_integrable_map,
         conserved_passive_scalar_map, eos_, stale_depth,
         *(lazy_passive_list_.get_list()));

      // increment stale_depth since the inner values have been updated
      // but the outer values have not
      stale_depth+=reconstructor->delayed_staling_rate();
    }
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::check_mesh_and_ghost_size_(Block *block) const noexcept
{
  Field field = block->data()->field();

  // Check that the mesh size is appropriate
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(field.field_id("density"),&gx,&gy,&gz);
  ASSERT("EnzoMethodMHDVlct::compute",
	 "Active zones on each block must be >= ghost depth.",
	 nx >= gx && ny >= gy && nz >= gz);

  // Check that the (cell-centered) ghost depth is large enough
  // Face-centered ghost could in principle be 1 smaller
  int min_ghost_depth = (half_dt_recon_->total_staling_rate() +
			 full_dt_recon_->total_staling_rate());
  ASSERT1("EnzoMethodMHDVlct::compute", "ghost depth must be at least %d.",
	  min_ghost_depth, std::min(nx, std::min(ny, nz)) >= min_ghost_depth);
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute_specific_passive_scalars_
(const str_vec_t &passive_list, EFlt3DArray& density,
 EnzoEFltArrayMap& conserved_passive_scalar_map,
 EnzoEFltArrayMap& specific_passive_scalar_map, int stale_depth) const noexcept
{
  int mz = density.shape(0);
  int my = density.shape(1);
  int mx = density.shape(2);

  for (const std::string& key : passive_list){
    EFlt3DArray cur_conserved = conserved_passive_scalar_map.at(key);
    EFlt3DArray out_specific = specific_passive_scalar_map.at(key);

    for (int iz = stale_depth; iz < mz - stale_depth; iz++) {
      for (int iy = stale_depth; iy < my - stale_depth; iy++) {
        for (int ix = stale_depth; ix < mx - stale_depth; ix++) {
          out_specific(iz,iy,ix) = cur_conserved(iz,iy,ix)/density(iz,iy,ix);
        }
      }
    }

  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute_flux_
(int dim, double cur_dt, enzo_float cell_width,
 EnzoEFltArrayMap &reconstructable_map,
 EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
 EFlt3DArray &pressure_l, EFlt3DArray &pressure_r,
 EnzoEFltArrayMap &flux_map, EnzoEFltArrayMap &dUcons_map,
 EFlt3DArray *interface_velocity_arr_ptr, EnzoReconstructor &reconstructor,
 EnzoBfieldMethod *bfield_method, int stale_depth,
 const str_vec_t& passive_list) const noexcept
{

  // purely for the purposes of making the caluclation more explicit, we define
  // the following aliases for priml_map and primr_map
  EnzoEFltArrayMap &reconstructable_l = priml_map;
  EnzoEFltArrayMap &reconstructable_r = primr_map;
  EnzoEFltArrayMap &integrable_l = priml_map;
  EnzoEFltArrayMap &integrable_r = primr_map;

  // First, reconstruct the left and right interface values
  reconstructor.reconstruct_interface(reconstructable_map,
                                      reconstructable_l, reconstructable_r,
				      dim, eos_, stale_depth, passive_list);

  // We temporarily increment the stale_depth for the rest of this calculation
  // here. We can't fully increment otherwise it will screw up the
  // reconstruction along other dimensions
  int cur_stale_depth = stale_depth + reconstructor.immediate_staling_rate();

  // Need to set the component of reconstructed B-field along dim, equal to
  // the corresponding longitudinal component of the B-field tracked at cell
  // interfaces (should potentially be handled internally by reconstructor)
  if (bfield_method != nullptr) {
    bfield_method->correct_reconstructed_bfield(reconstructable_l,
                                                reconstructable_r,
                                                dim, cur_stale_depth);
  }

  // Calculate integrable values on left and right faces:
  eos_->integrable_from_reconstructable(reconstructable_l, integrable_l,
					cur_stale_depth, passive_list);
  eos_->integrable_from_reconstructable(reconstructable_r, integrable_r,
					cur_stale_depth, passive_list);

  // Calculate pressure on left and right faces:
  eos_->pressure_from_reconstructable(reconstructable_l, pressure_l,
                                      cur_stale_depth);
  eos_->pressure_from_reconstructable(reconstructable_r, pressure_r,
                                      cur_stale_depth);

  // Next, compute the fluxes
  riemann_solver_->solve(integrable_l, integrable_r, pressure_l, pressure_r,
                         flux_map, dim, eos_, cur_stale_depth, passive_list,
                         interface_velocity_arr_ptr);

  // Accumulate the change in integrable quantities from these flux_map in
  // dUcons_map
  integrable_updater_->accumulate_flux_component(dim, cur_dt, cell_width,
                                                 flux_map, dUcons_map,
                                                 cur_stale_depth, passive_list);

  // if using dual energy formalism, compute the dual energy formalism
  if (eos_->uses_dual_energy_formalism()){
    EnzoSourceInternalEnergy eint_src;
    eint_src.calculate_source(dim, cur_dt, cell_width, reconstructable_map,
                              dUcons_map, *interface_velocity_arr_ptr, eos_,
                              cur_stale_depth);
  }

  // Finally, have the record the upwind direction (for handling CT)
  if (bfield_method != nullptr){
    bfield_method->identify_upwind(flux_map, dim, cur_stale_depth);
  }
}

//----------------------------------------------------------------------

static void add_temporary_arrays_to_map_
(EnzoEFltArrayMap &map, std::array<int,3> &shape,
 const std::vector<std::string>* const nonpassive_names,
 const str_vec_t* const passive_lists)
{

  if (nonpassive_names != nullptr){
    for (const std::string& name : (*nonpassive_names)){
      map[name] = EFlt3DArray(shape[0], shape[1], shape[2]);
    }
  }

  if (passive_lists != nullptr){
    for (const std::string& key : (*passive_lists)){
      map[key] = EFlt3DArray(shape[0], shape[1], shape[2]);
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::setup_arrays_
(Block *block, EnzoEFltArrayMap &primitive_map,
 EnzoEFltArrayMap &temp_primitive_map,
 EnzoEFltArrayMap &conserved_passive_scalar_map,
 EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
 EFlt3DArray &pressure_l, EFlt3DArray &pressure_r,
 EnzoEFltArrayMap &xflux_map, EnzoEFltArrayMap &yflux_map,
 EnzoEFltArrayMap &zflux_map, EnzoEFltArrayMap &dUcons_map) noexcept
{

  // allocate stuff! Make sure to do it in a way such that we don't have to
  // separately deallocate it!

  // To assist with setting up arrays, let's create a list of ALL primitive
  // keys (including passive scalars) and all integrable keys. These are the
  // same thing except the latter excludes quantities (like pressure) that we
  // don't compute pressure for.
  std::vector<std::string> combined_key_list = unique_combination_
    (integrable_field_list_, reconstructable_field_list_);

  // First, setup conserved_passive_scalar_map
  conserved_passive_scalar_map = conserved_passive_scalar_map_(block);

  // Next, setup nonpassive components of primitive_map
  primitive_map = nonpassive_primitive_map_(block);
  std::array<int,3> shape = {primitive_map.at("density").shape(0),
                             primitive_map.at("density").shape(1),
                             primitive_map.at("density").shape(2)};
  add_temporary_arrays_to_map_(primitive_map, shape, nullptr,
                               (lazy_passive_list_.get_list()).get());

  // Then, setup temp_primitive_map
  add_temporary_arrays_to_map_(temp_primitive_map, shape, &combined_key_list,
                               (lazy_passive_list_.get_list()).get());

  // Prepare arrays to hold fluxes (it should include groups for all actively
  // and passively advected quantities)
  EnzoEFltArrayMap* flux_maps[3] = {&zflux_map, &yflux_map, &xflux_map};
  for (std::size_t i = 0; i < 3; i++){
    std::array<int,3> cur_shape = shape; // makes a deep copy
    cur_shape[i] -= 1;
    add_temporary_arrays_to_map_(*(flux_maps[i]), cur_shape,
                                 &combined_key_list,
                                 (lazy_passive_list_.get_list()).get());
  }

  // Prepare fields used to accumulate all changes to the actively advected and
  // passively advected quantities. If CT is in use, dUcons_group should not
  // have storage for magnetic fields since CT independently updates magnetic
  // fields (this exclusion is implicitly handled integrable_updater_)
  std::vector<std::string> tmp = integrable_updater_->integrable_keys();
  add_temporary_arrays_to_map_(dUcons_map, shape, &tmp,
                               (lazy_passive_list_.get_list()).get());

  // Prepare temporary fields for priml and primr
  // As necessary, we pretend that these are centered along:
  //   - z and have shape (mz-1,  my,  mx)
  //   - y and have shape (  mz,my-1,  mx)
  //   - x and have shape (  mz,  my,mx-1)
  add_temporary_arrays_to_map_(priml_map, shape, &combined_key_list,
                               (lazy_passive_list_.get_list()).get());
  add_temporary_arrays_to_map_(primr_map, shape, &combined_key_list,
                               (lazy_passive_list_.get_list()).get());

  // If there are pressure entries in priml_map and primr_map (depends on the
  // EOS), set pressure_l and pressure_name_r equal to
  // those field names. Otherwise, reserve/allocate left/right pressure fields
  if (priml_map.contains("pressure")) {
    pressure_l = priml_map.at("pressure");
    pressure_r = primr_map.at("pressure");
  } else {
    pressure_l = EFlt3DArray(shape[0], shape[1], shape[2]);
    pressure_r = EFlt3DArray(shape[0], shape[1], shape[2]);
  }

}

//----------------------------------------------------------------------

double EnzoMethodMHDVlct::timestep ( Block * block ) const throw()
{
  // analogous to ppm timestep calulation, probably want to require that cfast
  // is no smaller than some tiny positive number.

  // Construct a map holding the field data for each of the (non-passive)
  // primitive quantities.
  EnzoEFltArrayMap primitive_map = nonpassive_primitive_map_(block);

  if (eos_->uses_dual_energy_formalism()){
    // synchronize eint and etot.
    // This is only strictly necessary after problem initialization and when
    // there is an inflow boundary condition
    eos_->apply_floor_to_energy_and_sync(primitive_map, 0);
  }

  // Compute thermal pressure (this presently requires that "pressure" is a
  // permanent field)
  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray pressure = array_factory.from_name("pressure");
  EnzoEFltArrayMap conserved_passive_scalar_map =
      conserved_passive_scalar_map_(block);
  eos_->pressure_from_integrable(primitive_map, pressure,
				 conserved_passive_scalar_map, 0);

  // Now load other necessary quantities
  enzo_float gamma = eos_->get_gamma();
  EFlt3DArray density = primitive_map.at("density");
  EFlt3DArray velocity_x = primitive_map.at("velocity_x");
  EFlt3DArray velocity_y = primitive_map.at("velocity_y");
  EFlt3DArray velocity_z = primitive_map.at("velocity_z");

  const bool mhd = (mhd_choice_ != bfield_choice::no_bfield);
  EFlt3DArray bfieldc_x, bfieldc_y, bfieldc_z;
  if (mhd) {
    bfieldc_x = primitive_map.at("bfield_x");
    bfieldc_y = primitive_map.at("bfield_y");
    bfieldc_z = primitive_map.at("bfield_z");
  }

  // widths of cells
  EnzoBlock * enzo_block = enzo::block(block);
  double dx = enzo_block->CellWidth[0];
  double dy = enzo_block->CellWidth[1];
  double dz = enzo_block->CellWidth[2];

  // initialize
  double dtBaryons = ENZO_HUGE_VAL;

  // timestep is the minimum of 0.5 * dr_i/(abs(v_i)+cfast) for all dimensions.
  // dr_i and v_i are the the width of the cell and velocity along dimension i.
  // cfast = fast magnetosonic speed (Convention is to use max value:
  // cfast = (va^2+cs^2)

  for (int iz=0; iz<density.shape(0); iz++) {
    for (int iy=0; iy<density.shape(1); iy++) {
      for (int ix=0; ix<density.shape(2); ix++) {
	enzo_float bmag_sq = 0.0;
        if (mhd){
          bmag_sq = (bfieldc_x(iz,iy,ix) * bfieldc_x(iz,iy,ix) +
		     bfieldc_y(iz,iy,ix) * bfieldc_y(iz,iy,ix) +
		     bfieldc_z(iz,iy,ix) * bfieldc_z(iz,iy,ix));
        }
	// Using "Rationalized" Gaussian units (where the magnetic permeability
	// is mu=1 and pressure = B^2/2)
	// To convert B to normal Gaussian units, multiply by sqrt(4*pi)
	enzo_float inv_dens= 1./density(iz,iy,ix);
	double cfast = (double) std::sqrt(gamma * pressure(iz,iy,ix) *
					  inv_dens + bmag_sq * inv_dens);

	dtBaryons = std::min(dtBaryons,
			     dx/(std::fabs((double) velocity_x(iz,iy,ix)) +
				 cfast));
	dtBaryons = std::min(dtBaryons,
			     dy/(std::fabs((double) velocity_y(iz,iy,ix)) +
				 cfast));
	dtBaryons = std::min(dtBaryons,
			     dz/(std::fabs((double) velocity_z(iz,iy,ix))
				 + cfast));
      }
    }
  }
  // Multiply resulting dt by CourantSafetyNumber (for extra safety!).
  // This should be less than 0.5 for standard algorithm
  dtBaryons *= courant_;

  return dtBaryons;
}
