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
				      double gamma, double density_floor,
				      double pressure_floor)
  : Method()
{
  // Initialize equation of state (check the validity of quantity floors)
  EnzoEquationOfState::check_floor(density_floor);
  EnzoEquationOfState::check_floor(pressure_floor);
  eos_ = new EnzoEOSIdeal(gamma, density_floor, pressure_floor);


  // determine integrable and reconstructable quantities (and passive scalars)
  determine_quantities_(eos_, integrable_group_names_,
			reconstructable_group_names_, passive_group_names_);

  // setup the groupings now
  setup_groupings_(integrable_group_names_, reconstructable_group_names_,
		   passive_group_names_);

  // "pressure" is only used to compute the timestep
  FieldDescr * field_descr = cello::field_descr();
  ASSERT("EnzoMethodMHDVlct", "\"pressure\" must be a permanent field",
	 field_descr->is_field("pressure"));

  // setup primitive_group_, and bfieldi_group_
  // (also checks that the integrable fields of primitive_group_ and all the
  // fields of bfieldi_group_ exist and are permanent)
  setup_groupings_(integrable_group_names_, reconstructable_group_names_,
		   passive_group_names_);

  // Initialize the default Refresh object - eventually may want to adjust
  // number of ghost zones based on reconstructor choice.
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
  			     enzo_sync_id_method_vlct);
  refresh(ir)->add_all_fields();

  // Should probably only refesh what we need
  // (using primitive_group_ and bfieldi_group_)
  //FieldDescr * field_descr = cello::field_descr();
  //refresh(ir)->add_field(field_descr->field_id("density"));
  //refresh(ir)->add_field(field_descr->field_id("velocity_x"));
  //refresh(ir)->add_field(field_descr->field_id("velocity_y"));
  //refresh(ir)->add_field(field_descr->field_id("velocity_z"));
  //refresh(ir)->add_field(field_descr->field_id("pressure"));
  //refresh(ir)->add_field(field_descr->field_id("bfield_x"));
  //refresh(ir)->add_field(field_descr->field_id("bfield_y"));
  //refresh(ir)->add_field(field_descr->field_id("bfield_z"));
  //refresh(ir)->add_field(field_descr->field_id("bfieldi_x"));
  //refresh(ir)->add_field(field_descr->field_id("bfieldi_y"));
  //refresh(ir)->add_field(field_descr->field_id("bfieldi_z"));

  // Initialize the component objects
  half_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (reconstructable_group_names_, passive_group_names_, half_recon_name);
  full_dt_recon_ = EnzoReconstructor::construct_reconstructor
    (reconstructable_group_names_, passive_group_names_, full_recon_name);
  riemann_solver_ = EnzoRiemann::construct_riemann
    (integrable_group_names_,      passive_group_names_, rsolver);
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::determine_quantities_
(EnzoEquationOfState *eos, std::vector<std::string> &integrable_quantities,
 std::vector<std::string> &reconstructable_quantities,
 std::vector<std::string> &passive_groups)
{
  if (enzo::config()->method_grackle_use_grackle){
    // This is in no way ready to use Grackle. Other things need to be dealt
    // with first

    // Need to make sure that all the passively advected Grackle fields are set
    // up. If the grackle method were to be called after this method, then some
    // passive fields might get setup that should be tracked by this method
    ERROR("EnzoMethodMHDVlct::determine_quantities_",
	  "Not presently equipped to handle grackle");
  }

  if (eos->uses_dual_energy_formalism()){
    ERROR("EnzoMethodMHDVlct::determine_quantities_",
	  "Not presently equipped to handle dual energy formalism.");
  }

  if (eos->is_barotropic()){
    ERROR("EnzoMethodMHDVlct::determine_quantities_",
	  "Not presently equipped to handle barotropic equations of state.");
  }
  
  std::string common[3] {"density", "velocity", "bfield"};
  for (int i = 0; i<3; i++){
    integrable_quantities.push_back(common[i]);
    reconstructable_quantities.push_back(common[i]);
  }

  // add specific total energy to integrable quantities
  integrable_quantities.push_back("total_energy");
  // add specific internal energy to reconstructable quantities
  reconstructable_quantities.push_back("internal_energy");

  // At a later date come back and deal with passive scalars

}

//----------------------------------------------------------------------

void add_passive_groups_(Grouping &grouping, const std::string field_prefix)
{
  // Helper function that adds entries of passively advected groups (specified
  // in the configuration file) an existing Grouping
  EnzoCenteredFieldRegistry registry;
  std::vector<std::string> group_names = registry.passive_scalar_group_names();
  Grouping* ref_grouping = cello::field_descr()->groups();
  for (unsigned int i=0;i<group_names.size();i++){
    std::string group_name = group_names[i];
    int num_fields = ref_grouping->size(group_name);
    for (int j=0;j<num_fields;j++){
      // Determine field_name
      std::string field_name = field_prefix + ref_grouping->item(group_name,j);
      grouping.add(field_name, group_name);
    }
  }
}

//----------------------------------------------------------------------

// Returns the unique members of a combination of 2 vectors
std::vector<std::string> unique_combination_(std::vector<std::string> &a,
					     std::vector<std::string> &b)
{
  std::vector<std::string> out;
  std::copy(a.begin(), a.end(), out.begin());
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

  // Add passive scalars to primitive_group_ (These are temporary specific
  // versions of the permanent fields)
  // This should be double checked
  //add_passive_groups_(*primitive_group_, "prim_specific_");
  ASSERT("EnzoMethodMHDVlct::setup_groupings_",
	 "Not quite expecting passive scalars", passive_groups.size() == 0);

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
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute ( Block * block) throw()
{
  if (block->is_leaf()) {
    // Check that the mesh size and ghost depths are appropriate
    check_mesh_and_ghost_size_(block);
    
    // declaring Grouping that track temporary fields used for scratch space
    // Groupings of conserved fields and primitive fields are tracked as members
    // conserved_group_ and primitive_group_

    // left and right reconstructed primitives
    Grouping priml_group, primr_group;

    // Names of temporary fields used to store the pressure computed from the
    // reconstructed left and right primitives
    std::string pressure_name_l, pressure_name_r;

    // flux fields
    Grouping xflux_group, yflux_group, zflux_group;

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
			  xflux_group, yflux_group, zflux_group,
			  efield_group, center_efield_name, weight_group,
			  temp_primitive_group, temp_bfieldi_group);

    // allocate constrained transport object
    EnzoConstrainedTransport ct = EnzoConstrainedTransport();

    double dt = block->dt();

    // convert the passive scalars from conserved form to specific form
    // (outside the integrator, they are treated like conserved densities)
    compute_specific_passive_scalars_(block, passive_group_names_,
				      *primitive_group_);

    // stale_depth indicates the number of field entries from the outermost
    // field value that the region including "stale" values (need to be
    // refreshed) extends over.
    int stale_depth = 0;

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
      }

      // Compute the reconstructable quantities from the integrable quantites
      // (note that primitve_groups_ actually holds all of the relevant
      // primitives since there is large overlap).
      // For an adiabatic gas, this nominally computes the specific internal
      // energy
      eos_->reconstructable_from_integrable(block, *cur_integrable_group,
					    *cur_reconstructable_group,
					    stale_depth);
      
      // Compute flux along each dimension
      compute_flux_(block, 0, *cur_reconstructable_group, *cur_bfieldi_group,
		    priml_group, primr_group, pressure_name_l, pressure_name_r,
		    xflux_group, weight_group, *reconstructor, stale_depth);
      compute_flux_(block, 1, *cur_reconstructable_group, *cur_bfieldi_group,
		    priml_group, primr_group, pressure_name_l, pressure_name_r,
		    yflux_group, weight_group, *reconstructor, stale_depth);
      compute_flux_(block, 2, *cur_reconstructable_group, *cur_bfieldi_group,
		    priml_group, primr_group, pressure_name_l, pressure_name_r,
		    zflux_group, weight_group, *reconstructor, stale_depth);
      // increment the stale_depth
      stale_depth+=reconstructor->immediate_staling_rate();

      // Compute the edge-centered Electric fields (each time, it uses the
      // initial cell-centered integrable quantities from the start of the
      // time-step
      compute_efields_(block, *primitive_group_,
		       xflux_group, yflux_group, zflux_group,
		       center_efield_name, efield_group,
		       weight_group, ct, stale_depth);

      // Add source terms? - this can happen after adding the flux divergence

      // but then the code placing a floor on the total energy must be moved.
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

      // Update quantities - add flux divergence
      // (this needs to happen updating the cell-centered B-field so that the
      //  floor on the energy can be checked)
      update_quantities_(block, *primitive_group_,
			 xflux_group, yflux_group, zflux_group,
			 *out_integrable_group, cur_dt, stale_depth);

      // increment stale_depth since the inner values have been updated
      // but the outer values have not
      stale_depth+=reconstructor->delayed_staling_rate();
    }

    // convert the passive scalars from specific form back to conserved form
    // (outside the integrator, they are treated like conserved densities)
    compute_conserved_passive_scalars_(block, passive_group_names_,
				       *primitive_group_, stale_depth);

    // Deallocate Temporary Fields
    deallocate_temp_fields_(block, priml_group, primr_group,
			    pressure_name_l, pressure_name_r,
			    xflux_group, yflux_group, zflux_group,
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


void EnzoMethodMHDVlct::compute_flux_(Block *block, int dim,
				      Grouping &reconstructable_group,
				      Grouping &cur_bfieldi_group,
				      Grouping &priml_group,
				      Grouping &primr_group,
				      std::string pressure_name_l,
				      std::string pressure_name_r,
				      Grouping &flux_group,
				      Grouping &weight_group,
				      EnzoReconstructor &reconstructor,
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
					*integrable_group_l, stale_depth, dim);
  eos_->integrable_from_reconstructable(block, *reconstructable_group_r,
					*integrable_group_r, stale_depth, dim);

  // Calculate pressure on left and right faces:
  eos_->pressure_from_reconstructable(block, *reconstructable_group_l,
				      pressure_name_l, stale_depth, dim);
  eos_->pressure_from_reconstructable(block, *reconstructable_group_r,
				      pressure_name_r, stale_depth, dim);

  // Next, compute the fluxes
  riemann_solver_->solve(block, *integrable_group_l, *integrable_group_r,
			 pressure_name_l, pressure_name_r, flux_group,
			 dim, eos_, cur_stale_depth);

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

void EnzoMethodMHDVlct::compute_efields_(Block *block,
					 Grouping &initial_integrable_group,
					 Grouping &xflux_group,
					 Grouping &yflux_group,
					 Grouping &zflux_group,
					 std::string center_efield_name,
					 Grouping &efield_group,
					 Grouping &weight_group,
					 EnzoConstrainedTransport &ct,
					 int stale_depth)
{
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // Maybe the following should be handled internally by ct?
  for (int i = 0; i < 3; i++){
    ct.compute_center_efield (block, i, center_efield_name,
			      initial_integrable_group, stale_depth);
    Grouping *jflux_group;
    Grouping *kflux_group;
    if (i == 0){
      jflux_group = &yflux_group;
      kflux_group = &zflux_group;
    } else if (i==1){
      jflux_group = &zflux_group;
      kflux_group = &xflux_group;
    } else {
      jflux_group = &xflux_group;
      kflux_group = &yflux_group;
    }

    ct.compute_edge_efield (block, i, center_efield_name, efield_group,
			    *jflux_group, *kflux_group, *primitive_group_,
			    weight_group, stale_depth);
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::update_quantities_(Block *block,
					   Grouping &initial_integrable_group,
					   Grouping &xflux_group,
					   Grouping &yflux_group,
					   Grouping &zflux_group,
					   Grouping &out_integrable_group,
					   double dt, int stale_depth)
{
  // Going to factor this out into a separate object
  std::vector<std::string> integrable_groups = {"density", "velocity",
						"total_energy", "bfield"};
  // Don't want to add fluxes to bfield since we use constrained transport
  std::vector<std::string> flagged_groups = {"bfield"};

  EnzoAdvectionFieldLUT lut;
  int conserved_start, conserved_stop;
  int specific_start, specific_stop;
  int other_start, other_stop;
  int nfields;

  EnzoCenteredFieldRegistry registry;
  lut = registry.prepare_advection_lut(integrable_groups,
				       conserved_start, conserved_stop,
				       specific_start, specific_stop,
				       other_start, other_stop, nfields,
				       flagged_groups);

  ASSERT("EnzoMethodMHDVlct::update_quantities_",
	 "not equipped to handle fields classified as other",
	 other_start == other_stop);

  EFlt3DArray *cur_prim_arrays, *out_prim_arrays;
  cur_prim_arrays = registry.load_array_of_fields(block, lut, nfields,
						  initial_integrable_group,
						  stale_depth);
  out_prim_arrays = registry.load_array_of_fields(block, lut, nfields,
						  out_integrable_group,
						  stale_depth);
  
  EFlt3DArray *xflux_arrays, *yflux_arrays, *zflux_arrays;
  xflux_arrays = registry.load_array_of_fields(block, lut, nfields,
					       xflux_group, stale_depth);
  yflux_arrays = registry.load_array_of_fields(block, lut, nfields,
					       yflux_group, stale_depth);
  zflux_arrays = registry.load_array_of_fields(block, lut, nfields,
					       zflux_group, stale_depth);

  enzo_float *cur_prim, *dU;
  cur_prim = new enzo_float[nfields];
  dU = new enzo_float[nfields];
  
  // For now, not having density floor affect momentum or total energy density
  enzo_float density_floor = eos_->get_density_floor();

  // cell-centered grid dimensions
  int mz = cur_prim_arrays[lut.density].shape(0);
  int my = cur_prim_arrays[lut.density].shape(1);
  int mx = cur_prim_arrays[lut.density].shape(2);

  // widths of cells
  EnzoBlock * enzo_block = enzo::block(block);
  enzo_float dtdx = dt/enzo_block->CellWidth[0];
  enzo_float dtdy = dt/enzo_block->CellWidth[1];
  enzo_float dtdz = dt/enzo_block->CellWidth[2];

  for (int iz=1; iz<mz-1; iz++) {
    for (int iy=1; iy<my-1; iy++) {
      for (int ix=1; ix<mx-1; ix++) {

	// load in the fields
	for (int field_ind=0; field_ind<nfields; field_ind++){
	  cur_prim[field_ind] = cur_prim_arrays[field_ind](iz,iy,ix);
	  dU[field_ind] = (- dtdx * (xflux_arrays[field_ind](iz,iy,ix)
				    - xflux_arrays[field_ind](iz,iy,ix-1))
			   - dtdy * (yflux_arrays[field_ind](iz,iy,ix)
				     - yflux_arrays[field_ind](iz,iy-1,ix))
			   - dtdz * (zflux_arrays[field_ind](iz,iy,ix)
				     - zflux_arrays[field_ind](iz-1,iy,ix)));
	}

	// get the initial density
	enzo_float initial_density = cur_prim[lut.density];

	// now update the integrable primitives that are conserved
	for (int i = conserved_start; i < conserved_stop; i++){
	  out_prim_arrays[i](iz,iy,ix) = cur_prim[i] + dU[i];
	}

	enzo_float new_density = out_prim_arrays[lut.density](iz,iy,ix);
	// possibly place a floor on new_density.
	new_density = EnzoEquationOfState::apply_floor(new_density,
						       density_floor);
	out_prim_arrays[lut.density](iz,iy,ix) = new_density;

	// update the specific primitives
	for (int i = specific_start; i < specific_stop; i++){
	  out_prim_arrays[i](iz,iy,ix) = (cur_prim[i] * initial_density
					  + dU[i]) /new_density;
	}

	// handle dual energy formalism ...

	// We are assuming that there aren't any quantities in the other
	// category (i.e. not classified as conserved/specific)

      }
    }
  }

  // apply floor to energy
  eos_->apply_floor_to_total_energy(block, out_integrable_group, stale_depth+1);

  // add fluxes to specific passive scalars (this could nominally be done in
  // the main update loop (if we were

  delete[] cur_prim;         delete[] dU;
  delete[] cur_prim_arrays;  delete[] out_prim_arrays;
  delete[] xflux_arrays;     delete[] yflux_arrays;     delete[] zflux_arrays;
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

// Helper function used to prepare groupings of temporary fields
// 
// In short, for every group-field combination in ref_grouping, where the group
// is also listed in group_names, an analogous group-field combination is
// added to grouping. The names of groups in ref_grouping and grouping are the
// same while the temporary fields in grouping are given by field_prefix + the
// name of the analogous field in ref_grouping. Additionally, all temporary
// fields added to grouping are allocated with centering given by cx, cy, cz.
void prep_temp_field_grouping_(Field &field, Grouping &ref_grouping,
			       std::vector<std::string> &group_names,
			       Grouping &grouping, std::string field_prefix,
			       int cx, int cy, int cz)
{
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

void EnzoMethodMHDVlct::allocate_temp_fields_(Block *block,
					      Grouping &priml_group,
					      Grouping &primr_group,
					      std::string &pressure_name_l,
					      std::string &pressure_name_r,
					      Grouping &xflux_group,
					      Grouping &yflux_group,
					      Grouping &zflux_group,
					      Grouping &efield_group,
					      std::string &center_efield_name,
					      Grouping &weight_group,
					      Grouping &temp_primitive_group,
					      Grouping &temp_bfieldi_group)
{
  
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();
  FieldDescr * field_descr = field.field_descr();

  ASSERT("EnzoMethodMHDVlct::allocate_temp_fields_",
	 "Need to allocate temporary fields to hold specific passive scalars.",
	 passive_group_names_.size() == 0);
  
  // First, reserve/allocate temporary fields for any relevant reconstructable
  // primitives that are not also integrable and have not yet been allocated
  for (unsigned int i=0;i<reconstructable_group_names_.size();i++){
    std::string group_name = reconstructable_group_names_[i];
    int num_fields = primitive_group_->size(group_name);

    for (int j=0;j<num_fields;j++){
      // Determine field_name
      std::string field_name = primitive_group_->item(group_name,j);

      if (field.is_field(field_name) &&
	  field.is_permanent(field.field_id(field_name))){
	continue;
      }

      ASSERT1("EnzoMethodMHDVlct::allocate_temp_fields_",
	      "%s is integrable, it should be a permanent field.",
	      field_name.c_str(),
	      (std::find(integrable_group_names_.begin(),
			 integrable_group_names_.end(), group_name) !=
	       integrable_group_names_.end()));

      // allocate temporary field
      int ir_ = field_descr->insert_temporary(field_name);
      field.allocate_temporary(ir_);
    }
  }

  // Prepare the temporary primitive fields (used to store values at half dt)
  std::vector<std::string> prim_group_names;
  prim_group_names = unique_combination_(integrable_group_names_,
					 reconstructable_group_names_);
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names,
			    temp_primitive_group, "temp_",0,0,0);

  // Prepare temporary flux fields
  prep_temp_field_grouping_(field, *primitive_group_, integrable_group_names_,
			    xflux_group, "xflux_",-1,0,0);
  prep_temp_field_grouping_(field, *primitive_group_, integrable_group_names_,
			    yflux_group, "yflux_",0,-1,0);
  prep_temp_field_grouping_(field, *primitive_group_, integrable_group_names_,
			    zflux_group, "zflux_",0,0,-1);

  // Prepare temporary fields for priml and primr
  // As necessary, we pretend that these are centered along:
  //   - z and have shape (mz-1,  my,  mx)
  //   - y and have shape (  mz,my-1,  mx)
  //   - x and have shape (  mz,  my,mx-1)
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names,
			    priml_group, "left_",0,0,0);
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names,
			    primr_group, "right_",0,0,0);

  // reserve/allocate left/right pressure fields
  pressure_name_l = "left_single_pressure";
  pressure_name_r = "right_single_pressure";
  prep_reused_temp_field_(field, pressure_name_l, 0, 0, 0);
  prep_reused_temp_field_(field, pressure_name_r, 0, 0, 0);

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

  // I'm commenting out the following. With the introduction of stale depth
  // the following should no longer be necessary

  /*
  // Initialize Values (that would not otherwise be initialized) in order to
  // avoid NaN/Inf during calculations.

  // initialize values of the temporary conserved group (to avoid NaNs during
  // conversion to primitives in ghost zone - these propogate to flux calc) 
  outer_ghost_to_floor_(block, temp_conserved_group, *eos_, cons_group_names_);
  
  // initialize exterior face values of the temporary interface B-field
  // (to avoid NaNs while computing cell-centered B-fields in ghost zone)
  std::vector<std::string> bfieldi_group_names{"bfield"};
  outer_ghost_to_floor_(block,temp_bfieldi_group, *eos_, bfieldi_group_names);

  // initialize values of reconstructed primitives
  initialize_recon_prim_to_floor_(block, priml_group, *eos_, prim_group_names_);
  initialize_recon_prim_to_floor_(block, primr_group, *eos_, prim_group_names_);
  */
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

void EnzoMethodMHDVlct::deallocate_temp_fields_(Block *block,
						Grouping &priml_group,
						Grouping &primr_group,
						std::string pressure_name_l,
						std::string pressure_name_r,
						Grouping &xflux_group,
						Grouping &yflux_group,
						Grouping &zflux_group,
						Grouping &efield_group,
						std::string center_efield_name,
						Grouping &weight_group,
						Grouping &temp_primitive_group,
						Grouping &temp_bfieldi_group)
{

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  ASSERT("EnzoMethodMHDVlct::deallocate_temp_fields_",
	 ("Need to deallocate temporary fields holding specific passive "
	  "scalars."),
	 passive_group_names_.size() == 0);
  
  // Prepare the temporary primitive fields (used to store values at half dt)
  std::vector<std::string> prim_group_names;
  prim_group_names = unique_combination_(integrable_group_names_,
					 reconstructable_group_names_);

  // deallocate cell-centered primitive groups
  deallocate_grouping_fields_(field, prim_group_names, *primitive_group_);
  deallocate_grouping_fields_(field, prim_group_names, temp_primitive_group);
  deallocate_grouping_fields_(field, prim_group_names, priml_group);
  deallocate_grouping_fields_(field, prim_group_names, primr_group);
  deallocate_grouping_fields_(field, integrable_group_names_, xflux_group);
  deallocate_grouping_fields_(field, integrable_group_names_, yflux_group);
  deallocate_grouping_fields_(field, integrable_group_names_, zflux_group);

  // deallocate left/right pressure fields
  field.deallocate_temporary(field.field_id(pressure_name_l));
  field.deallocate_temporary(field.field_id(pressure_name_r));
  

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

void EnzoMethodMHDVlct::initialize_recon_prim_to_floor_
(Block *block, Grouping &grouping, EnzoEquationOfState &eos,
 std::vector<std::string> &group_names)
{
  int start = -1;
  int stop = -1;
  Field field = block->data()->field();

  for (unsigned int group_ind=0;group_ind<group_names.size(); group_ind++){
    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = grouping.size(group_name);

    // Handle possibility of having a density/pressure floor
    enzo_float prim_floor = 0.;
    if (group_name == "density"){
      prim_floor = eos.get_density_floor();
    } else if (group_name == "pressure"){
      prim_floor = eos.get_pressure_floor();
    }

    for (int field_ind=0; field_ind<num_fields; field_ind++){
      std::string field_name = grouping.item(group_name,field_ind);
      if (start == -1){
	// Determine the indices to set to floor. Only need to initialize
	// values at the end of the allocated blocks of memory.
	int id = field.field_id(field_name);
	int mx, my, mz;
	field.dimensions(id,&mx,&my,&mz);
	stop = mz*my*mx;
	start = std::min(std::min(mz*my*(mx-1), mz*(my-1)*mx),(mz-1)*my*mx);
      }
      enzo_float* data = (enzo_float *) field.values(field_name);
      for (int i=start; i<stop; i++) {
	data[i] = prim_floor;
      }
    }
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::outer_ghost_to_floor_
(Block *block, Grouping &grouping, EnzoEquationOfState &eos,
 std::vector<std::string> &group_names){

  EnzoFieldArrayFactory array_factory(block);
  
  for (unsigned int group_ind=0;group_ind<group_names.size();group_ind++){
    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = grouping.size(group_name);


    // Handle possibility of having a density/pressure floor
    enzo_float floor_val = 0.;
    if (group_name == "density"){
      floor_val = eos.get_density_floor();
    } else if (group_name == "pressure"){
      floor_val = eos.get_pressure_floor();
    } else if (group_name == "total_energy"){
      // We actually primarily expect to be doing this for conserved quantities
      floor_val = eos.get_pressure_floor() / (eos.get_gamma() - 1.);
    } 
     // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      // load in the quantities
      EFlt3DArray array;
      array = array_factory.from_grouping(grouping, group_name, field_ind);
      int mz = array.shape(0);
      int my = array.shape(1);
      int mx = array.shape(2);
 
      // set values for the left edges
      array.subarray(CSlice(0,1), CSlice(0,my), CSlice(0,mx)) = floor_val;
      array.subarray(CSlice(0,mz), CSlice(0,1), CSlice(0,mx)) = floor_val;
      array.subarray(CSlice(0,mz), CSlice(0,my), CSlice(0,1)) = floor_val;

      // set values for the right edges
      array.subarray(CSlice(-1,mz), CSlice(0,my), CSlice(0,mx)) = floor_val;
      array.subarray(CSlice(0,mz), CSlice(-1,my), CSlice(0,mx)) = floor_val;
      array.subarray(CSlice(0,mz), CSlice(0,my), CSlice(-1,mx)) = floor_val;
    }
  }
}

//----------------------------------------------------------------------

double EnzoMethodMHDVlct::timestep ( Block * block ) const throw()
{
  // analogous to ppm timestep calulation, probably want to require that cfast
  // is no smaller than some tiny positive number.

  // Compute the pressure (assumes that "pressure" is a permanent field)
  FieldDescr * field_descr = cello::field_descr();
  eos_->pressure_from_integrable(block, *primitive_group_, "pressure",
				 *(field_descr->groups()), false, 0);
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

  // The start index of the active region
  int ix_start = enzo_block->GridStartIndex[0];
  int iy_start = enzo_block->GridStartIndex[1];
  int iz_start = enzo_block->GridStartIndex[2];

  // The end index (final valid index) of the active region
  int ix_end = enzo_block->GridEndIndex[0];
  int iy_end = enzo_block->GridEndIndex[1];
  int iz_end = enzo_block->GridEndIndex[2];

  // widths of cells
  double dx = enzo_block->CellWidth[0];
  double dy = enzo_block->CellWidth[1];
  double dz = enzo_block->CellWidth[2];

  // initialize
  enzo_float dtBaryons = ENZO_HUGE_VAL;

  
  
  // timestep is the minimum of 0.5 * dr_i/(abs(v_i)+cfast) for all dimensions.
  // dr_i and v_i are the the width of the cell and velocity along dimension i.
  // cfast = fast magnetosonic speed (Convention is to use max value:
  // cfast = (va^2+cs^2)

  for (int iz=iz_start; iz<=iz_end; iz++) {
    for (int iy=iy_start; iy<=iy_end; iy++) {
      for (int ix=ix_start; ix<=ix_end; ix++) {
	enzo_float bmag_sq = (bfieldc_x(iz,iy,ix) * bfieldc_x(iz,iy,ix) +
			      bfieldc_y(iz,iy,ix) * bfieldc_y(iz,iy,ix) +
			      bfieldc_z(iz,iy,ix) * bfieldc_z(iz,iy,ix));

	// Using "Rationalized" Gaussian units (where the magnetic permeability
	// is mu=1 and pressure = B^2/2)
	// To convert B to normal Gaussian units, multiply by sqrt(4*pi)
	enzo_float inv_dens= 1./density(iz,iy,ix);
	enzo_float cfast = std::sqrt(gamma * pressure(iz,iy,ix) * inv_dens + 
				     bmag_sq * inv_dens);

	dtBaryons = std::min(dtBaryons,
			     dx/(std::fabs(velocity_x(iz,iy,ix)) + cfast));
	dtBaryons = std::min(dtBaryons,
			     dy/(std::fabs(velocity_y(iz,iy,ix)) + cfast));
	dtBaryons = std::min(dtBaryons,
			     dz/(std::fabs(velocity_z(iz,iy,ix)) + cfast));
      }
    }
  }
  /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
  dtBaryons *= courant_;

  ASSERT2("EnzoMethodMHDVlct::timestep",
	  "Invalid timestep, %g, was calculated for %s.",
	  (double)dtBaryons, block->name().c_str(),
	  ((dtBaryons>0) && (std::isfinite(dtBaryons))));
  fflush(stdout);

  return dtBaryons;
}
