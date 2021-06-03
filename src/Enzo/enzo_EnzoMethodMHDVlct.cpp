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

#ifdef CONFIG_USE_GRACKLE
  if (enzo::config()->method_grackle_use_grackle){
    // we can remove the following once EnzoMethodGrackle no longer requires
    // the internal_energy to be a permanent field
    ASSERT("EnzoMethodMHDVlct::determine_quantities_",
           ("Grackle cannot currently be used alongside this integrator "
            "unless the dual-energy formalism is in use"),
           eos_->uses_dual_energy_formalism());
  }
#endif /* CONFIG_USE_GRACKLE */

  // Determine whether magnetic fields are to be used
  mhd_choice_ = parse_bfield_choice_(mhd_choice);

  riemann_solver_ = EnzoRiemann::construct_riemann
    (rsolver, mhd_choice_ != bfield_choice::no_bfield,
     eos_->uses_dual_energy_formalism());

  // determine integration and primitive quantities
  std::vector<std::string> integration_quantities, primitive_quantities;
  integration_quantities = riemann_solver_->integration_quantities();
  primitive_quantities = riemann_solver_->primitive_quantities();

  // Initialize the remaining component objects
  half_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (primitive_quantities, half_recon_name, (enzo_float)theta_limiter);
  full_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (primitive_quantities, full_recon_name, (enzo_float)theta_limiter);
  
  integration_quan_updater_ =
    new EnzoIntegrationQuanUpdate(integration_quantities, true);

  // Determine the lists of fields that are required to hold the integration
  // quantities and primitives and ensure that they are defined
  build_field_l_(integration_quantities, integration_field_list_);
  build_field_l_(primitive_quantities, primitive_field_list_);

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

EnzoMethodMHDVlct::~EnzoMethodMHDVlct()
{
  delete eos_;
  delete half_dt_recon_;
  delete full_dt_recon_;
  delete riemann_solver_;
  delete integration_quan_updater_;
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
  p|integration_quan_updater_;
  p|mhd_choice_;
  p|bfield_method_;
  p|integration_field_list_;
  p|primitive_field_list_;
  p|lazy_passive_list_;
}

//----------------------------------------------------------------------

EnzoEFltArrayMap EnzoMethodMHDVlct::get_integration_map_
(Block * block,  const str_vec_t *passive_list) const noexcept
{
  EnzoEFltArrayMap map("integration");
  EnzoFieldArrayFactory array_factory(block,0);
  for (const std::string& field_name : integration_field_list_){
    ASSERT1("EnzoMethodMHDVlct::get_integration_map_",
            "EnzoEFltArrayMap can't hold more than one key called \"%s\"",
            field_name.c_str(), !map.contains(field_name));
    map[field_name] = array_factory.from_name(field_name);
  }

  if (passive_list != nullptr){
    for (const std::string& field_name : (*passive_list)){
      ASSERT1("EnzoMethodMHDVlct::get_integration_map_",
              "EnzoEFltArrayMap can't hold more than one key called \"%s\"",
              field_name.c_str(), !map.contains(field_name));
      map[field_name] = array_factory.from_name(field_name);
    }
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
    // integration quantities. Additionally, this also includes temporary
    // arrays used to hold the specific form of the passive scalar
    EnzoEFltArrayMap integration_map; // this will be overwritten

    // map used for storing integration values at the half time-step. This
    // includes key,array pairs for each entry in integration_map (there should
    // be no aliased fields shared between maps)
    EnzoEFltArrayMap temp_integration_map("temp_integration");

    // Map of arrays used to temporarily store the cell-centered primitive
    // quantities that are subsequently reconstructed. This includes arrays for
    // storing the specific form of each of the passively advected scalars.
    EnzoEFltArrayMap primitive_map("primitive");

    // map holding the arrays wrapping each fields corresponding to a passively
    // advected scalar. The scalar is in conserved form.
    EnzoEFltArrayMap conserved_passive_scalar_map; // this will be overwritten

    // holds left and right reconstructed primitives (scratch-space)
    EnzoEFltArrayMap priml_map("priml");
    EnzoEFltArrayMap primr_map("primr");

    // maps used to store fluxes (in the future, these will wrap FluxData
    // entries)
    EnzoEFltArrayMap xflux_map("xflux");
    EnzoEFltArrayMap yflux_map("yflux");
    EnzoEFltArrayMap zflux_map("zflux");

    // map of arrays  used to accumulate the changes to the conserved forms of
    // the integration quantities and passively advected scalars. In other
    // words, at the start of the (partial) timestep, the fields are all set to
    // zero and are used to accumulate the flux divergence and source terms. If
    // CT is used, it won't have space to store changes in the magnetic fields.
    EnzoEFltArrayMap dUcons_map("dUcons");

    setup_arrays_(block, integration_map, temp_integration_map, primitive_map,
                  priml_map, primr_map,  xflux_map, yflux_map, zflux_map,
                  dUcons_map);

    // Setup a pointer to an array that used to store interface velocity fields
    // from computed by the Riemann Solver (to use in the calculation of the
    // internal energy source term). If the dual energy formalism is not in
    // use, don't actually allocate the array and set the pointer to NULL.
    EFlt3DArray interface_velocity_arr, *interface_velocity_arr_ptr;
    if (eos_->uses_dual_energy_formalism()){
      EFlt3DArray density = integration_map.at("density");
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

    // repeat the following loop twice (for half time-step and full time-step)
    for (int i=0;i<2;i++){
      double cur_dt = (i == 0) ? dt/2. : dt;
      EnzoEFltArrayMap& cur_integration_map =
        (i == 0) ? integration_map      : temp_integration_map;
      EnzoEFltArrayMap& out_integration_map =
        (i == 0) ? temp_integration_map :      integration_map;

      EnzoReconstructor *reconstructor;

      if (i == 0){
        reconstructor = half_dt_recon_;
      } else {
        reconstructor = full_dt_recon_;
      }

      // set all elements of the arrays in dUcons_map to 0 (throughout the rest
      // of the current loop, flux divergence and source terms will be
      // accumulated in these arrays)
      integration_quan_updater_->clear_dUcons_map
        (dUcons_map, 0., *(lazy_passive_list_.get_list()));

      // Compute the primitive quantities from the integration quantites
      // This basically copies all quantities that are both and an integration
      // quantity and a primitive and converts the passsive scalars from
      // conserved-form to specific-form (i.e. from density to mass fraction).
      // For a non-barotropic gas, this also computes pressure
      eos_->primitive_from_integration(cur_integration_map, primitive_map,
                                       stale_depth,
                                       *(lazy_passive_list_.get_list()));

      // Compute flux along each dimension
      compute_flux_(0, cur_dt, cell_widths[0], primitive_map,
                    priml_map, primr_map, xflux_map, dUcons_map,
                    interface_velocity_arr_ptr, *reconstructor, bfield_method_,
                    stale_depth, *(lazy_passive_list_.get_list()));
      compute_flux_(1, cur_dt, cell_widths[1], primitive_map,
                    priml_map, primr_map, yflux_map, dUcons_map,
                    interface_velocity_arr_ptr, *reconstructor, bfield_method_,
                    stale_depth, *(lazy_passive_list_.get_list()));
      compute_flux_(2, cur_dt, cell_widths[2], primitive_map,
                    priml_map, primr_map, zflux_map, dUcons_map,
                    interface_velocity_arr_ptr, *reconstructor, bfield_method_,
                    stale_depth, *(lazy_passive_list_.get_list()));

      // increment the stale_depth
      stale_depth+=reconstructor->immediate_staling_rate();

      // This is where source terms should be computed (added to dUcons_group)

      // Update Bfields
      if (bfield_method_ != nullptr) {
        bfield_method_->update_all_bfield_components(cur_integration_map,
                                                     xflux_map, yflux_map,
                                                     zflux_map,
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
        (integration_map, dUcons_map, out_integration_map, eos_, stale_depth,
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

void EnzoMethodMHDVlct::compute_flux_
(const int dim, const double cur_dt, const enzo_float cell_width,
 EnzoEFltArrayMap &primitive_map,
 EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
 EnzoEFltArrayMap &flux_map, EnzoEFltArrayMap &dUcons_map,
 EFlt3DArray *interface_velocity_arr_ptr, EnzoReconstructor &reconstructor,
 EnzoBfieldMethod *bfield_method, const int stale_depth,
 const str_vec_t& passive_list) const noexcept
{

  // First, reconstruct the left and right interface values
  reconstructor.reconstruct_interface(primitive_map, priml_map, primr_map,
				      dim, eos_, stale_depth, passive_list);

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
  riemann_solver_->solve(priml_map, primr_map, flux_map, dim, eos_,
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
  if (eos_->uses_dual_energy_formalism()){
    EnzoSourceInternalEnergy eint_src;
    eint_src.calculate_source(dim, cur_dt, cell_width, primitive_map,
                              dUcons_map, *interface_velocity_arr_ptr, eos_,
                              cur_stale_depth);
  }

  // Finally, have bfield_method record the upwind direction (for handling CT)
  if (bfield_method != nullptr){
    bfield_method->identify_upwind(flux_map, dim, cur_stale_depth);
  }
}

//----------------------------------------------------------------------

static void add_temporary_arrays_to_map_
(EnzoEFltArrayMap &map, const std::array<int,3> &shape,
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
(Block *block, EnzoEFltArrayMap &integration_map,
 EnzoEFltArrayMap &temp_integration_map, EnzoEFltArrayMap &primitive_map,
 EnzoEFltArrayMap &priml_map, EnzoEFltArrayMap &primr_map,
 EnzoEFltArrayMap &xflux_map, EnzoEFltArrayMap &yflux_map,
 EnzoEFltArrayMap &zflux_map, EnzoEFltArrayMap &dUcons_map) noexcept
{

  // allocate stuff! Make sure to do it in a way such that we don't have to
  // separately deallocate it!

  // First, setup integration_map
  integration_map = get_integration_map_
    (block, (lazy_passive_list_.get_list()).get());

  const std::array<int,3> shape = {integration_map.at("density").shape(0),
				   integration_map.at("density").shape(1),
				   integration_map.at("density").shape(2)};

  // Next, setup temp_integration_map
  add_temporary_arrays_to_map_(temp_integration_map, shape,
                               &integration_field_list_,
                               (lazy_passive_list_.get_list()).get());

  // Prepare arrays to hold fluxes. It should include keys for all integration
  // quantities actively (including passively advected scalars)
  EnzoEFltArrayMap* flux_maps[3] = {&zflux_map, &yflux_map, &xflux_map};
  for (std::size_t i = 0; i < 3; i++){
    std::array<int,3> cur_shape = shape; // makes a deep copy
    cur_shape[i] -= 1;
    add_temporary_arrays_to_map_(*(flux_maps[i]), cur_shape,
                                 &integration_field_list_,
                                 (lazy_passive_list_.get_list()).get());
  }

  // Prepare fields used to accumulate all changes to the integration
  // quantities (including passively advected scalars). If CT is in use,
  // dUcons_group should not have storage for magnetic fields since CT
  // independently updates magnetic fields (this exclusion is implicitly
  // handled by integration_quan_updater_)
  std::vector<std::string> tmp = integration_quan_updater_->integration_keys();
  add_temporary_arrays_to_map_(dUcons_map, shape, &tmp,
                               (lazy_passive_list_.get_list()).get());

  // Setup primitive_map
  add_temporary_arrays_to_map_(primitive_map, shape,
                               &primitive_field_list_,
                               (lazy_passive_list_.get_list()).get());

  // Prepare maps for holding the left and right reconstructed primitives.
  // As necessary, we pretend that these are centered along:
  //   - z and have shape (mz-1,  my,  mx)
  //   - y and have shape (  mz,my-1,  mx)
  //   - x and have shape (  mz,  my,mx-1)
  add_temporary_arrays_to_map_(priml_map, shape, &primitive_field_list_,
                               (lazy_passive_list_.get_list()).get());
  add_temporary_arrays_to_map_(primr_map, shape, &primitive_field_list_,
                               (lazy_passive_list_.get_list()).get());
}

//----------------------------------------------------------------------

double EnzoMethodMHDVlct::timestep ( Block * block ) const throw()
{
  // analogous to ppm timestep calulation, probably want to require that cfast
  // is no smaller than some tiny positive number.

  // Constructs a map containing the field data for each integration quantity
  // This includes each passively advected scalar (as densities)
  EnzoEFltArrayMap integration_map = get_integration_map_
    (block, (lazy_passive_list_.get_list()).get());

  if (eos_->uses_dual_energy_formalism()){
    // synchronize eint and etot.
    // This is only strictly necessary after problem initialization and when
    // there is an inflow boundary condition
    eos_->apply_floor_to_energy_and_sync(integration_map, 0);
  }

  // Compute thermal pressure (this presently requires that "pressure" is a
  // permanent field)
  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray pressure = array_factory.from_name("pressure");
  eos_->pressure_from_integration(integration_map, pressure, 0);

  // Now load other necessary quantities
  const enzo_float gamma = eos_->get_gamma();
  EFlt3DArray density = integration_map.at("density");
  EFlt3DArray velocity_x = integration_map.at("velocity_x");
  EFlt3DArray velocity_y = integration_map.at("velocity_y");
  EFlt3DArray velocity_z = integration_map.at("velocity_z");

  const bool mhd = (mhd_choice_ != bfield_choice::no_bfield);
  EFlt3DArray bfieldc_x, bfieldc_y, bfieldc_z;
  if (mhd) {
    bfieldc_x = integration_map.at("bfield_x");
    bfieldc_y = integration_map.at("bfield_y");
    bfieldc_z = integration_map.at("bfield_z");
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
