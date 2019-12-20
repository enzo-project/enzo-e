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

EnzoMethodMHDVlct::EnzoMethodMHDVlct (std::string rsolver,
				      std::string half_recon_name,
				      std::string full_recon_name,
				      double gamma, double theta_limiter,
				      double density_floor,
				      double pressure_floor,
				      bool dual_energy_formalism,
				      double dual_energy_formalism_eta)
  : Method()
{
  // Initialize equation of state (check the validity of quantity floors)
  EnzoEquationOfState::check_floor(density_floor);
  EnzoEquationOfState::check_floor(pressure_floor);
  eos_ = new EnzoEOSIdeal(gamma, density_floor, pressure_floor,
			  dual_energy_formalism, dual_energy_formalism_eta);


  // determine integrable and reconstructable quantities (and passive scalars)
  determine_quantities_(eos_, integrable_group_names_,
			reconstructable_group_names_, passive_group_names_);

  FieldDescr * field_descr = cello::field_descr();

  // "pressure" is only used to compute the timestep
  ASSERT("EnzoMethodMHDVlct", "\"pressure\" must be a permanent field",
	 field_descr->is_field("pressure"));


  // setup primitive_group_, and bfieldi_group_
  // (also checks that the integrable fields of primitive_group_ and all the
  // fields of bfieldi_group_ exist and are permanent)
  setup_groupings_(integrable_group_names_, reconstructable_group_names_,
		   passive_group_names_);


  // Initialize the default Refresh object - May want to adjust
  // number of ghost zones based on reconstructor choice.
  const int ir = add_refresh(4, 0, neighbor_leaf, sync_barrier,
  			     enzo_sync_id_method_vlct);
  // Add cell-centered fields holding actively advected integrable quantities
  add_group_fields_to_refresh_(ir, *primitive_group_, integrable_group_names_);
  // Add cell-centered fields holding the conserved form of the passively
  // advected quantities
  add_group_fields_to_refresh_(ir, *(field_descr->groups()),
			       passive_group_names_);

  refresh(ir)->add_field(field_descr->field_id("bfieldi_x"));
  refresh(ir)->add_field(field_descr->field_id("bfieldi_y"));
  refresh(ir)->add_field(field_descr->field_id("bfieldi_z"));

  // Initialize the remaining component objects
  half_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (reconstructable_group_names_, passive_group_names_, half_recon_name,
     (enzo_float)theta_limiter);
  full_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (reconstructable_group_names_, passive_group_names_, full_recon_name,
     (enzo_float)theta_limiter);
  riemann_solver_ = EnzoRiemann::construct_riemann
    (integrable_group_names_,      passive_group_names_, rsolver);
  integrable_updater_ = new EnzoIntegrableUpdate(integrable_group_names_,
						 true, passive_group_names_);
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::determine_quantities_
(EnzoEquationOfState *eos, std::vector<std::string> &integrable_quantities,
 std::vector<std::string> &reconstructable_quantities,
 std::vector<std::string> &passive_groups)
{
  if (enzo::config()->method_grackle_use_grackle){
    // make sure all the required fields are defined so that the group of
    // "colour" fields is accurate (needed for identifying passive scalars)
    EnzoMethodGrackle::define_required_grackle_fields();
    // Not quite ready to support a variable gamma
  }


  std::string common[3] {"density", "velocity", "bfield"};
  for (int i = 0; i<3; i++){
    integrable_quantities.push_back(common[i]);
    reconstructable_quantities.push_back(common[i]);
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

  // Now to setup the list of passively advected group names

  // retrieve list of all possible group names of passively advected scalars
  EnzoCenteredFieldRegistry registry;
  std::vector<std::string> all_passive_group_names;
  all_passive_group_names = registry.passive_scalar_group_names();

  // We'll now iterate through all known possible groups names that could
  // contain passive scalars and check to see if any fields were labelled as
  // being a part of these groups in the input file (if so, then add the group
  // name to the list of passively advected group names)
  Grouping* reference_grouping = cello::field_descr()->groups();
  for (std::size_t i=0; i < all_passive_group_names.size(); i++){
    std::string group_name = all_passive_group_names[i];

    if (reference_grouping->size(group_name) > 0){
      passive_groups.push_back(group_name);
    }
  }

}

//----------------------------------------------------------------------

void add_passive_groups_(Grouping &grouping, const std::string field_prefix,
			 const std::vector<std::string> group_names)
{
  // Helper function that adds entries of passively advected groups (specified
  // in the configuration file) to an existing Grouping object

  // This includes the groups specified in the parameter file
  Grouping* reference_grouping = cello::field_descr()->groups();

  for (unsigned int i=0;i<group_names.size();i++){

    std::string group_name = group_names[i];
    int num_fields = reference_grouping->size(group_name);

    for (int j=0; j < num_fields; j++){
      // Determine field_name
      std::string field_name = (field_prefix +
				reference_grouping->item(group_name,j));

      // add the field to the grouping
      grouping.add(field_name, group_name);
    }
  }
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

void EnzoMethodMHDVlct::setup_groupings_
(std::vector<std::string> &integrable_groups,
 std::vector<std::string> &reconstructable_groups,
 std::vector<std::string> &passive_groups)
{

  FieldDescr * field_descr = cello::field_descr();

  // first come up with a vector group names that represents the union of
  // integrable_groups and reconstructable_groups
  std::vector<std::string> groups;
  groups = unique_combination_(integrable_groups,reconstructable_groups);

  // now setup primitive_group_ using the names in groups
  EnzoCenteredFieldRegistry registry;
  primitive_group_ = registry.build_grouping(groups, "");

  // We should check that all the fields in integrable groups are real
  // permenant fields
  for (std::size_t i = 0; i<integrable_groups.size(); i++){
    std::string group_name = integrable_groups[i];
    int num_fields = primitive_group_->size(group_name);

    for (int j = 0; j<num_fields; j++){
      std::string field_name = primitive_group_->item(group_name,j);

      ASSERT1("EnzoMethodMHDVlct::setup_groupings_",
	      "\"%s\" must be the name of a permanent field",
	      field_name.c_str(), field_descr->is_field(field_name));
    }
  }

  // Add temporary fields to primitive_group_ to hold specific forms (mass
  // fraction) of passively advected scalars. These get are used to hold the
  // values after they are converted from conserved form (density) which are
  // tracked in permanent fields
  add_passive_groups_(*primitive_group_, "specific_passive_", passive_groups);

  bfieldi_group_ = new Grouping;
  bfieldi_group_->add("bfieldi_x", "bfield");
  bfieldi_group_->add("bfieldi_y", "bfield");
  bfieldi_group_->add("bfieldi_z", "bfield");

  ASSERT("EnzoMethodMHDVlct::setup_groupings_",
	 ("There must be face-centered permanent fields called \"bfieldi_x\","
	  "\"bfield_y\" and \"bfield_z\"."),
	 (field_descr->is_field("bfieldi_x") &&
	  field_descr->is_field("bfieldi_y") &&
	  field_descr->is_field("bfieldi_z")));
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::add_group_fields_to_refresh_
(const int ir, Grouping &grouping, std::vector<std::string> group_names)
{
  FieldDescr * field_descr = cello::field_descr();

  for (unsigned int i=0;i<group_names.size();i++){

    std::string group_name = group_names[i];
    int num_fields = grouping.size(group_name);

    for (int j=0; j < num_fields; j++){
      // Determine field_name
      std::string field_name = grouping.item(group_name,j);
      ASSERT1("EnzoMethodMHDVlct::add_group_fields_to_refresh_",
	      "%s must be a permanent field",field_name.c_str(),
	      field_descr->is_permanent(field_name))

      // add the field to the refresh object
      refresh(ir)->add_field(field_descr->field_id(field_name));
    }
  }
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::~EnzoMethodMHDVlct()
{
  delete primitive_group_;
  delete bfieldi_group_;

  delete eos_;
  delete half_dt_recon_;
  delete full_dt_recon_;
  delete riemann_solver_;
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p|eos_;
  p|integrable_group_names_;
  p|reconstructable_group_names_;
  p|passive_group_names_;

  if (p.isUnpacking()){
    // a bug prevents us from pupping the groupings
    // instead, we just set them up again
    setup_groupings_(integrable_group_names_, reconstructable_group_names_,
		     passive_group_names_);
  }

  p|half_dt_recon_;
  p|full_dt_recon_;
  p|riemann_solver_;
  p|integrable_updater_;
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute ( Block * block) throw()
{
  if (block->is_leaf()) {
    // Check that the mesh size and ghost depths are appropriate
    check_mesh_and_ghost_size_(block);
    
    // declaring Groupings that track temporary fields used for scratch space
    // left and right reconstructed primitives
    Grouping priml_group, primr_group;

    // Names of temporary fields used to store the pressure computed from the
    // reconstructed left and right primitives
    // Note: in the case of adiabatic fluids, pressure is a reconstructable
    //       quantity and fields are included for it in priml_group and
    //       primr_group. In that case, those field names are also asigned to
    //       pressure_name_l and pressure_name_r (this is not strictly,
    //       necessary - EnzoEOSIdeal is smart enough to copy values if the
    //       field names don't match)
    std::string pressure_name_l, pressure_name_r;

    // Name of fields use to store interface velocity field from the Riemann
    // Solver to compute internal energy source term. It won't be allocated
    // unless the dual energy formalism is in use
    std::string interface_velocity_name;

    // flux fields
    Grouping xflux_group, yflux_group, zflux_group;

    // fields used to accumulate the changes to the conserved forms of the
    // integrable quantities and passively advected scalars (i.e. at the start
    // of the (partial) timestep, the fields are all set to zero and are used
    // to accumulate the flux divergence and source terms). If CT is used, it
    // won't have space to store changes in the magnetic fields.
    Grouping dUcons_group;

    // edge-centered electric fields
    Grouping efield_group;

    // Name of the temporary field to store the central E-field (can reuse it
    // to store different components of the field at different times)
    std::string center_efield_name;

    // face-centered weight fields - Track the x/y/z upwind direction
    Grouping weight_group;

    // temp primitive group for storing values at the half time-step
    Grouping temp_primitive_group;

    // temp interface b-fileds
    Grouping temp_bfieldi_group;

    // allocate the temporary fields (as necessary) and fill the field groupings
    allocate_temp_fields_(block, priml_group, primr_group,
			  pressure_name_l, pressure_name_r,
			  interface_velocity_name,
			  xflux_group, yflux_group, zflux_group, dUcons_group,
			  efield_group, center_efield_name, weight_group,
			  temp_primitive_group, temp_bfieldi_group);

    // conserved-form of passively advected scalars (the initial values are
    // taken from here and converted to specific form. Intermediate and final
    // values get written to here)
    Grouping *conserved_passive_scalars = cello::field_descr()->groups();

    // allocate constrained transport object
    EnzoConstrainedTransport ct = EnzoConstrainedTransport();

    double dt = block->dt();

    // stale_depth indicates the number of field entries from the outermost
    // field value that the region including "stale" values (need to be
    // refreshed) extends over.
    int stale_depth = 0;

    // convert the passive scalars from conserved form to specific form
    // (outside the integrator, they are treated like conserved densities)
    compute_specific_passive_scalars_(block, passive_group_names_,
				      *conserved_passive_scalars,
				      *primitive_group_, stale_depth);

    // repeat the following loop twice (for half time-step and full time-step)

    for (int i=0;i<2;i++){
      double cur_dt;
      // For the purposes of making the calculation procedure slightly more
      // more explicit, we distinguish between cur_integrable_group and
      // cur_reconstructable_group. Due to the high level of overlap between
      // these, they are simply aliases of the same underlying Grouping that
      // holds groups for both of them
      Grouping* cur_reconstructable_group;
      Grouping* cur_integrable_group;
      Grouping* out_integrable_group;

      Grouping* cur_bfieldi_group;
      Grouping* out_bfieldi_group;
      EnzoReconstructor *reconstructor;

      if (i == 0){
	cur_dt = dt/2.;
	out_integrable_group = &temp_primitive_group;
	cur_integrable_group = primitive_group_;
	cur_reconstructable_group = cur_integrable_group;
	out_bfieldi_group = &temp_bfieldi_group;
	cur_bfieldi_group = bfieldi_group_;
	reconstructor = half_dt_recon_;
      } else {
	cur_dt = dt;
	out_integrable_group = primitive_group_;
	cur_integrable_group = &temp_primitive_group;
	cur_reconstructable_group = cur_integrable_group;
	out_bfieldi_group = bfieldi_group_;
	cur_bfieldi_group = &temp_bfieldi_group;
	reconstructor = full_dt_recon_;

	// After the fluxes were added to the passive scalar in the first half
	// timestep, the values were stored in conserved form in the fields
	// held by conserved_passive_scalars.
	// Need to convert them to specific form
	compute_specific_passive_scalars_(block, passive_group_names_,
					  *conserved_passive_scalars,
					  *cur_integrable_group, stale_depth);
      }

      // set the all values of fields in dUcons_group to 0 (throughout the rest
      // of the current loop, flux divergence and source terms will be
      // accumulated by these fields)
      integrable_updater_->clear_dUcons_group(block, dUcons_group, 0.);

      // Compute the reconstructable quantities from the integrable quantites
      // Although cur_integrable_group holds the passive scalars in integrable
      // form, the conserved form of the values is required in case Grackle is
      // being used.
      //
      // Note: cur_integrable_group and cur_reconstructable_group refer to the
      // same underlying grouping since there is such a large degree of overlap
      // between reconstructable and integrable quantities
      //
      // For a barotropic gas, the following nominally does nothing
      // For a non-barotropic gas, the following nominally computes pressure
      eos_->reconstructable_from_integrable(block, *cur_integrable_group,
					    *cur_reconstructable_group,
					    *conserved_passive_scalars,
					    stale_depth);

      // Compute flux along each dimension
      compute_flux_(block, 0, cur_dt, *cur_reconstructable_group,
		    *cur_bfieldi_group, priml_group, primr_group,
		    pressure_name_l, pressure_name_r, xflux_group,
		    weight_group, dUcons_group, interface_velocity_name,
		    *reconstructor, stale_depth);
      compute_flux_(block, 1, cur_dt, *cur_reconstructable_group,
		    *cur_bfieldi_group, priml_group, primr_group,
		    pressure_name_l, pressure_name_r, yflux_group,
		    weight_group, dUcons_group, interface_velocity_name,
		    *reconstructor, stale_depth);
      compute_flux_(block, 2, cur_dt, *cur_reconstructable_group,
		    *cur_bfieldi_group, priml_group, primr_group,
		    pressure_name_l, pressure_name_r, zflux_group,
		    weight_group, dUcons_group, interface_velocity_name,
		    *reconstructor, stale_depth);
      // increment the stale_depth
      stale_depth+=reconstructor->immediate_staling_rate();

      // This is where source terms should be computed (added to dUcons_group)

      // Update Bfields (this could probably all be consolidated)
      // First, compute the edge-centered Electric fields (each time, it uses
      // the current integrable quantities
      ct.compute_all_edge_efields(block, *cur_integrable_group,
				  xflux_group, yflux_group, zflux_group,
				  center_efield_name, efield_group,
				  weight_group, stale_depth);

      // Update longitudinal B-field (add source terms of constrained transport)
      for (int dim = 0; dim<3; dim++){
	ct.update_bfield(block, dim, efield_group, *bfieldi_group_,
			 *out_bfieldi_group, cur_dt, stale_depth);
      }

      // Finally, update cell-centered B-field 
      for (int dim = 0; dim<3; dim++){
	ct.compute_center_bfield(block, dim, *out_integrable_group,
				 *out_bfieldi_group, stale_depth);
      }

      // Update quantities - (includes flux divergence and source terms) 
      // (this needs to happen after updating the cell-centered B-field so that
      // the pressure floor can be applied to the total energy (and if
      // necessary the total energy can be synchronized with internal energy)
      //
      // Note: updated passive scalars are NOT saved in out_integrable_group in
      //     specific form. Instead they are saved inconserved_passive_scalars
      //     in conserved form.
      integrable_updater_->update_quantities(block, *primitive_group_,
					     dUcons_group,
					     *out_integrable_group,
					     *conserved_passive_scalars,
					     eos_, stale_depth);

      // increment stale_depth since the inner values have been updated
      // but the outer values have not
      stale_depth+=reconstructor->delayed_staling_rate();
    }

    // Deallocate Temporary Fields
    deallocate_temp_fields_(block, priml_group, primr_group,
			    pressure_name_l, pressure_name_r,
			    interface_velocity_name, xflux_group,
			    yflux_group, zflux_group, dUcons_group,
			    efield_group, center_efield_name,
			    weight_group, temp_primitive_group,
			    temp_bfieldi_group);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::check_mesh_and_ghost_size_(Block *block)
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
  (Block *block, const std::vector<std::string> passive_groups,
   Grouping &conserved_passive_scalars,
   Grouping &primitive_group, int stale_depth)
{
  EnzoFieldArrayFactory array_factory(block, stale_depth);

  EFlt3DArray density = array_factory.from_grouping(primitive_group, "density",
						    0);

  std::vector<std::string> group_names = this->passive_group_names_;

  for (std::size_t group_ind = 0; group_ind<group_names.size(); group_ind++){

    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = conserved_passive_scalars.size(group_name);

    ASSERT("EnzoMethodMHDVlct::compute_specific_passive_scalars_",
	   "There shouldn't be any passive groups without any fields.",
	   num_fields != 0);

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){
      // load values
      EFlt3DArray cur_conserved, out_specific;
      cur_conserved = array_factory.from_grouping(conserved_passive_scalars,
						  group_name, field_ind);
      out_specific = array_factory.from_grouping(primitive_group, group_name,
						 field_ind);

      for (int iz=0; iz<density.shape(0); iz++) {
	for (int iy=0; iy<density.shape(1); iy++) {
	  for (int ix=0; ix<density.shape(2); ix++) {
	    out_specific(iz,iy,ix) = cur_conserved(iz,iy,ix)/density(iz,iy,ix);
	  }
	}
      }

    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute_flux_
(Block *block, int dim, double cur_dt, Grouping &reconstructable_group,
 Grouping &cur_bfieldi_group, Grouping &priml_group, Grouping &primr_group,
 std::string pressure_name_l, std::string pressure_name_r,
 Grouping &flux_group, Grouping &weight_group, Grouping &dUcons_group,
 std::string interface_velocity_name, EnzoReconstructor &reconstructor,
 int stale_depth)
{
  // purely for the purposes of making the caluclation more explicit, we define
  // the following aliases for priml_group/primr_group
  Grouping *integrable_group_l, *integrable_group_r;
  Grouping *reconstructable_group_l, *reconstructable_group_r;
  integrable_group_l = &priml_group;  reconstructable_group_l = &priml_group;
  integrable_group_r = &primr_group;  reconstructable_group_r = &primr_group;
  
  // First, reconstruct the left and right interface values
  reconstructor.reconstruct_interface(block, reconstructable_group,
				      *reconstructable_group_l,
				      *reconstructable_group_r,
				      dim, eos_, stale_depth);

  // We temporarily increment the stale_depth for the rest of this calculation
  // here. We can't fully increment otherwise it will screw up the
  // reconstruction along other dimensions
  int cur_stale_depth = stale_depth + reconstructor.immediate_staling_rate();

  // Need to set the component of reconstructed B-field along dim, equal to
  // the corresponding longitudinal component of the B-field tracked at cell
  // interfaces (should probably be handled internally by reconstructor)
  EnzoFieldArrayFactory array_factory(block, cur_stale_depth);
  EFlt3DArray bfield, l_bfield,r_bfield;
  bfield = array_factory.interior_bfieldi(cur_bfieldi_group, dim);
  l_bfield = array_factory.reconstructed_field(*reconstructable_group_l,
					       "bfield", dim, dim);
  r_bfield = array_factory.reconstructed_field(*reconstructable_group_r,
					       "bfield", dim, dim);

  // All 3 array objects are the same shape
  for (int iz=0; iz<bfield.shape(0); iz++) {
    for (int iy=0; iy<bfield.shape(1); iy++) {
      for (int ix=0; ix<bfield.shape(2); ix++) {
	l_bfield(iz,iy,ix) = bfield(iz,iy,ix);
	r_bfield(iz,iy,ix) = bfield(iz,iy,ix);
      }
    }
  }
  // Equivalently: l_bfield.subarray() = bfield;
  //               r_bfield.subarray() = bfield;

  // Calculate integrable values on left and right faces:
  eos_->integrable_from_reconstructable(block, *reconstructable_group_l,
					*integrable_group_l,
					cur_stale_depth, dim);
  eos_->integrable_from_reconstructable(block, *reconstructable_group_r,
					*integrable_group_r,
					cur_stale_depth, dim);

  // Calculate pressure on left and right faces:
  eos_->pressure_from_reconstructable(block, *reconstructable_group_l,
				      pressure_name_l,
				      cur_stale_depth, dim);
  eos_->pressure_from_reconstructable(block, *reconstructable_group_r,
				      pressure_name_r,
				      cur_stale_depth, dim);

  // Next, compute the fluxes
  riemann_solver_->solve(block, *integrable_group_l, *integrable_group_r,
			 pressure_name_l, pressure_name_r, flux_group,
			 dim, eos_, cur_stale_depth, interface_velocity_name);

  // Accumulate the change in integrable quantities from these fluxes in
  // dUcons_group
  integrable_updater_->accumulate_flux_component(block, dim, cur_dt,
						 flux_group, dUcons_group,
						 cur_stale_depth);

  // if using dual energy formalism, compute the dual energy formalism
  if (eos_->uses_dual_energy_formalism()){
    EnzoSourceInternalEnergy eint_src;
    eint_src.calculate_source(block, cur_dt, reconstructable_group,
			      dUcons_group, interface_velocity_name, dim, eos_,
			      cur_stale_depth);
  }

  // Finally, need to handle weights
  //  - Currently, weight is set to 1.0 if upwind is in positive direction of
  //    the current dimension, 0 if upwind is in the negative direction of the
  //    current dimension, or 0.5 if there is no upwind direction
  //  - At present, the weights are unnecessary (the same information is
  //    encoded in density flux to figure out this information). However, this
  //    functionallity is implemented in case we decide to adopt the weighting
  //    scheme from Athena++, which requires knowledge of the reconstructed
  //    densities.

  EFlt3DArray density_flux, weight_field;
  density_flux = array_factory.from_grouping(flux_group, "density", 0);
  weight_field = array_factory.from_grouping(weight_group, "weight", dim);

  // Iteration limits compatible with both 2D and 3D grids
  for (int iz=0; iz<density_flux.shape(0); iz++) {
    for (int iy=0; iy<density_flux.shape(1); iy++) {
      for (int ix=0; ix<density_flux.shape(2); ix++) {
	// density flux is the face-centered density times the face-centered
	// velocity along dim

	if ( density_flux(iz,iy,ix) > 0){
	  weight_field(iz,iy,ix) = 1.0;
	} else if ( density_flux(iz,iy,ix) < 0){
	  weight_field(iz,iy,ix) = 0.0;
	} else {
	  weight_field(iz,iy,ix) = 0.5;
	}
      }
    }
  }

}

//----------------------------------------------------------------------

// Helper function to prepare temporary fields that are reused each cycle
void prep_reused_temp_field_(Field &field, std::string field_name,
			     int cx, int cy, int cz)
{
  FieldDescr * field_descr = field.field_descr();
  int id;
  if (field.is_field(field_name)){
    // the temporary field was added in a previous time-step
    id = field.field_id(field_name);
  } else {
    // Reserve a temporary field
    id = field_descr->insert_temporary(field_name);
    field_descr->set_centering(id,cx,cy,cz);
  }

  // allocate temporary field
  // this last check allows us to make temporary field "permanent" for debugging
  // purposes
  if (!field.is_permanent(id)){
    field.allocate_temporary(id);
  }
}
//----------------------------------------------------------------------

void allocate_temp_group_fields_(Field &field, Grouping &grouping,
				 std::vector<std::string> &group_names,
				 int cx, int cy, int cz)
{
  // allocate the temporary fields in grouping
  for (std::size_t i = 0; i < group_names.size(); i++){
    std::string group_name = group_names[i];
    int num_fields = grouping.size(group_name);
    for (int j = 0; j < num_fields; j++){
      prep_reused_temp_field_(field, grouping.item(group_name,j), cx, cy, cz);
    }
  }
}

//----------------------------------------------------------------------

// Helper function used to prepare groupings of temporary fields
// 
// In short, for every group-field combination in ref_grouping, where the group
// is also listed in group_names, an analogous group-field combination is
// added to grouping. The names of groups in ref_grouping and grouping are the
// same while the temporary fields in grouping are given by field_prefix + the
// name of the analogous field in ref_grouping. Additionally, all temporary
// fields added to grouping are allocated with centering given by cx, cy, cz.
void prep_temp_field_grouping_(Field &field, Grouping &ref_grouping,
			       const std::vector<std::string> &group_names,
			       Grouping &grouping, std::string field_prefix,
			       int cx, int cy, int cz)
{
  // we can totally come up with an iterator to factor out the following nested
  // for loop (it appears a lot!)
  for (unsigned int i=0;i<group_names.size();i++){
    std::string group_name = group_names[i];
    int num_fields = ref_grouping.size(group_name);
    if (num_fields == 0) {
      continue; // This group is not tracked by ref_grouping
    }

    for (int j=0;j<num_fields;j++){
      // Determine field_name
      std::string field_name = field_prefix + (ref_grouping.item(group_name,j));
      // reserve/allocate field
      prep_reused_temp_field_(field, field_name, cx, cy, cz);

      // add the temporary field to grouping
      grouping.add(field_name,group_name);

    }
  }
}

//----------------------------------------------------------------------

// Helper function used to prepare groupings of temporary fields. These
// groupings only contain 3 fields - corresponding to vector components
//  - These vectors are all expected to be face-centered or edge centered.
//    (with the direction of centering depending on the dimension of the
//     component)
//  - face-centered: [face_center = True] fields are all centered on the face
//    corresponding to the direction of the vector component (e.g. x component
//    of weights is centered on the faces along the x dimension)
//  - edge-centered: [face_center = False] fields are centered on the edge
//    corresponding to the direction not pointed to by the vector component
//    (e.g. z-component of edge E-field is centered on the x and y edges)
//  - If exterior_faces is true, then the field will include space for the
//    exterior faces. If it is false, then the field will not include space for
//    the exterior faces.
//  - the names of the temporary fields are given by field_prefix + x,
//    field_prefix + y and field_prefix + z

void prep_temp_vector_grouping_(Field &field, std::string group_name,
				Grouping &grouping, std::string field_prefix,
				bool face_center, bool exterior_faces)
{
  for (int i=0;i<3;i++){
    // prepare field name
    std::string field_name;
    int cx, cy, cz, delta;
    if (face_center) {
      cx = 0; cy = 0; cz = 0;
      delta = (exterior_faces) ? 1 : -1;
    } else {
      if (exterior_faces){
	cx = 1; cy = 1; cz = 1;
	delta = -1;
      } else {
	cx = -1; cy = -1; cz = -1;
	delta = 1;
      }
    }

    if (i == 0){
      cx += delta;
      field_name = field_prefix + "x";
    } else if (i == 1) {
      cy += delta;
      field_name = field_prefix + "y";
    } else {
      cz += delta;
      field_name = field_prefix + "z";
    }

    // reserve/allocate field
    prep_reused_temp_field_(field, field_name, cx, cy, cz);

    grouping.add(field_name,group_name);
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::allocate_temp_fields_
(Block *block, Grouping &priml_group, Grouping &primr_group,
 std::string &pressure_name_l, std::string &pressure_name_r,
 std::string &interface_velocity_name,
 Grouping &xflux_group, Grouping &yflux_group, Grouping &zflux_group,
 Grouping &dUcons_group, Grouping &efield_group,
 std::string &center_efield_name, Grouping &weight_group,
 Grouping &temp_primitive_group, Grouping &temp_bfieldi_group)
{
  
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();
  
  // First, reserve/allocate temporary fields for any relevant reconstructable
  // primitives that have not yet been allocated (i.e. they are not also
  // integrable primitives)
  allocate_temp_group_fields_(field, *primitive_group_,
			      reconstructable_group_names_, 0, 0, 0);
  // Next, reserve/allocate temporary fields to hold the specific form of all
  // passively advected primitives
  allocate_temp_group_fields_(field, *primitive_group_, passive_group_names_,
			      0, 0, 0);

  // Prepare temporary flux fields (it should include groups for all actively
  // and passively advected quantities)
  std::vector<std::string> combined_integrable_names;
  combined_integrable_names = unique_combination_(integrable_group_names_,
						  passive_group_names_);
  prep_temp_field_grouping_(field, *primitive_group_, combined_integrable_names,
			    xflux_group, "xflux_",-1,0,0);
  prep_temp_field_grouping_(field, *primitive_group_, combined_integrable_names,
			    yflux_group, "yflux_",0,-1,0);
  prep_temp_field_grouping_(field, *primitive_group_, combined_integrable_names,
			    zflux_group, "zflux_",0,0,-1);

  // Prepare fields used to accumulate all changes to the activel advected and
  // passively advected quantities. If CT is in use, dUcons_group should not
  // have storage for magnetic fields since CT independently updates magnetic
  // fields (this exlusion is implicitly handled integrable_updater_)
  prep_temp_field_grouping_
    (field, *primitive_group_,
     integrable_updater_->combined_integrable_groups(true),
     dUcons_group, "dUcons_",0,0,0);

  // Prepare the temporary primitive fields (used to store values at half dt)
  std::vector<std::string> prim_group_names;
  prim_group_names = unique_combination_(combined_integrable_names,
					 reconstructable_group_names_);
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names,
			    temp_primitive_group, "temp_",0,0,0);

  // Prepare temporary fields for priml and primr
  // As necessary, we pretend that these are centered along:
  //   - z and have shape (mz-1,  my,  mx)
  //   - y and have shape (  mz,my-1,  mx)
  //   - x and have shape (  mz,  my,mx-1)
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names,
			    priml_group, "left_",0,0,0);
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names,
			    primr_group, "right_",0,0,0);

  // If there are pressure fields in priml_group and primr_group (depends on 
  // the equation of state), set pressure_name_l and pressure_name_r equal to
  // those field names. Otherwise, reserve/allocate left/right pressure fields
  if (priml_group.size("pressure") != 0){
    pressure_name_l = priml_group.item("pressure",0);
    pressure_name_r = primr_group.item("pressure",0);
  } else {
    pressure_name_l = "left_single_pressure";
    pressure_name_r = "right_single_pressure";
    prep_reused_temp_field_(field, pressure_name_l, 0, 0, 0);
    prep_reused_temp_field_(field, pressure_name_r, 0, 0, 0);
  }

  // if the dual energy formalism is in use, allocate a fields to temporarily
  // store velocity normal to the cell interface at the cell interface. The
  // interface velocities are used to compute the internal energy source term.
  if (eos_->uses_dual_energy_formalism()){
    interface_velocity_name = "temp_cell_interface_velocity";
    prep_reused_temp_field_(field, interface_velocity_name, 0, 0, 0);
  } else {
    interface_velocity_name = "";
  }

  // reserve/allocate fields for weight fields
  // these are face-centered fields that store the upwind/downwind direction
  prep_temp_vector_grouping_(field, "weight", weight_group,
			     "temp_weight_", true, false);

  // reserve allocate temporary interface bfields fields (includes the exterior
  // faces of the grid)
  prep_temp_vector_grouping_(field, "bfield", temp_bfieldi_group,
			     "temp_bfieldi_",true, true);

  // allocate temporary efield fields
  // reserve/allocate fields for edge-centered electric fields
  prep_temp_vector_grouping_(field, "efield", efield_group,
			     "temp_efield_",false,false);

  // reserve/allocate cell-centered e-field
  center_efield_name = "center_efield";
  prep_reused_temp_field_(field, center_efield_name, 0, 0, 0);
}

//----------------------------------------------------------------------

// Helper function that deallocates the temporary fields of a grouping
void deallocate_grouping_fields_(Field &field,
				 std::vector<std::string> &group_names,
				 Grouping &grouping)
{
  for (unsigned int i=0;i<group_names.size();i++){
    std::string group_name = group_names[i];
    int num_fields = grouping.size(group_name);
    for (int j=0;j<num_fields;j++){
      // Determine field_name and id
      std::string field_name = (grouping.item(group_name,j));
      int id = field.field_id(field_name);
      // if it's a temporary field, deallocate it
      if (!field.is_permanent(id)) {
	field.deallocate_temporary(id);
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::deallocate_temp_fields_
(Block *block, Grouping &priml_group, Grouping &primr_group,
 std::string pressure_name_l, std::string pressure_name_r,
 std::string interface_velocity_name,
 Grouping &xflux_group, Grouping &yflux_group, Grouping &zflux_group,
 Grouping &dUcons_group, Grouping &efield_group, std::string center_efield_name,
 Grouping &weight_group, Grouping &temp_primitive_group,
 Grouping &temp_bfieldi_group)
{

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // Prepare the temporary primitive fields (used to store values at half dt)
  std::vector<std::string> combined_integrable_names, prim_group_names;
  combined_integrable_names = unique_combination_(integrable_group_names_,
						  passive_group_names_);
  prim_group_names = unique_combination_(combined_integrable_names,
					 reconstructable_group_names_);

  // deallocate cell-centered primitive groups
  deallocate_grouping_fields_(field, prim_group_names, *primitive_group_);
  deallocate_grouping_fields_(field, prim_group_names, temp_primitive_group);
  deallocate_grouping_fields_(field, prim_group_names, priml_group);
  deallocate_grouping_fields_(field, prim_group_names, primr_group);

  // deallocate face-centered flux fields
  deallocate_grouping_fields_(field, combined_integrable_names, xflux_group);
  deallocate_grouping_fields_(field, combined_integrable_names, yflux_group);
  deallocate_grouping_fields_(field, combined_integrable_names, zflux_group);

  // deallocate left/right pressure fields
  field.deallocate_temporary(field.field_id(pressure_name_l));
  field.deallocate_temporary(field.field_id(pressure_name_r));

  // deallocate fields used to hold interace velocities (if it was allocated)
  if (interface_velocity_name!=""){
    field.deallocate_temporary(field.field_id(interface_velocity_name));
  }

  // deallocate electric fields
  std::vector<std::string> efield_group_names{"efield"};
  deallocate_grouping_fields_(field, efield_group_names, efield_group);
  field.deallocate_temporary(field.field_id(center_efield_name));

  // deallocate weights
  std::vector<std::string> weight_group_names{"weight"};
  deallocate_grouping_fields_(field, weight_group_names, weight_group);

  // deallocate the temporary longitudinal bfields
  std::vector<std::string> bfieldi_group_names{"bfield"};
  deallocate_grouping_fields_(field, bfieldi_group_names, temp_bfieldi_group);
}

//----------------------------------------------------------------------

double EnzoMethodMHDVlct::timestep ( Block * block ) const throw()
{
  // analogous to ppm timestep calulation, probably want to require that cfast
  // is no smaller than some tiny positive number.

  // Compute the pressure (requires that "pressure" is a permanent field)
  FieldDescr * field_descr = cello::field_descr();

  if (eos_->uses_dual_energy_formalism()){
    // synchronize eint and etot.
    // This is only strictly necessary after problem initialization and when
    // there is an inflow boundary condition
    eos_->apply_floor_to_energy_and_sync(block, *primitive_group_, 0);
  }
  
  eos_->pressure_from_integrable(block, *primitive_group_, "pressure",
				 *(field_descr->groups()), 0);
  enzo_float gamma = eos_->get_gamma();

  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray density, velocity_x, velocity_y, velocity_z, pressure;
  EFlt3DArray bfieldc_x, bfieldc_y, bfieldc_z;
  density = array_factory.from_grouping(*primitive_group_, "density", 0);
  velocity_x = array_factory.from_grouping(*primitive_group_, "velocity", 0);
  velocity_y = array_factory.from_grouping(*primitive_group_, "velocity", 1);
  velocity_z = array_factory.from_grouping(*primitive_group_, "velocity", 2);
  bfieldc_x = array_factory.from_grouping(*primitive_group_, "bfield", 0);
  bfieldc_y = array_factory.from_grouping(*primitive_group_, "bfield", 1);
  bfieldc_z = array_factory.from_grouping(*primitive_group_, "bfield", 2);

  pressure = array_factory.from_name("pressure");

  // Get iteration limits
  // Like ppm and ppml, access active region info from enzo_block attributes
  EnzoBlock * enzo_block = enzo::block(block);

  // widths of cells
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
	enzo_float bmag_sq = (bfieldc_x(iz,iy,ix) * bfieldc_x(iz,iy,ix) +
			      bfieldc_y(iz,iy,ix) * bfieldc_y(iz,iy,ix) +
			      bfieldc_z(iz,iy,ix) * bfieldc_z(iz,iy,ix));

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

  ASSERT2("EnzoMethodMHDVlct::timestep",
	  "Invalid timestep, %g, was calculated for %s.",
	  (double)dtBaryons, block->name().c_str(),
	  ((dtBaryons>0) && (std::isfinite(dtBaryons))));
  fflush(stdout);

  return dtBaryons;
}
