
#include "cello.hpp"
#include "enzo.hpp"
#include "charm_enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodMHDVlct::EnzoMethodMHDVlct (std::string rsolver,
				      std::string half_recon_name,
				      std::string full_recon_name,
				      double gamma, double density_floor,
				      double pressure_floor)
  : Method()
{
  // Initialize the default Refresh object - eventually may want to adjust
  // number of ghost zones based on reconstructor choice.

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
  			       enzo_sync_id_method_vlct);

  // Not clear if all fields should be refreshed (like in PPML constructor)
  FieldDescr * field_descr = cello::field_descr();
  
  refresh(ir)->add_field(field_descr->field_id("density"));
  refresh(ir)->add_field(field_descr->field_id("velocity_x"));
  refresh(ir)->add_field(field_descr->field_id("velocity_y"));
  refresh(ir)->add_field(field_descr->field_id("velocity_z"));
  refresh(ir)->add_field(field_descr->field_id("pressure"));
  refresh(ir)->add_field(field_descr->field_id("bfield_x"));
  refresh(ir)->add_field(field_descr->field_id("bfield_y"));
  refresh(ir)->add_field(field_descr->field_id("bfield_z"));
  refresh(ir)->add_field(field_descr->field_id("bfieldi_x"));
  refresh(ir)->add_field(field_descr->field_id("bfieldi_y"));
  refresh(ir)->add_field(field_descr->field_id("bfieldi_z"));

  // Setup the EnzoFieldConditions struct
  EnzoFieldConditions cond;
  cond.hydro = true;
  cond.MHD = true;

  cond_ = cond;

  // setup primitive_group_, conserved_group_, bfieldi_group_ AND get group
  // names
  setup_groups_(cond);

  // Initialize the component objects
  eos_ = new EnzoEOSIdeal(gamma, density_floor, pressure_floor);
  half_dt_recon_ = EnzoReconstructor::construct_reconstructor(half_recon_name,
							      cond);
  full_dt_recon_ = EnzoReconstructor::construct_reconstructor(full_recon_name,
							      cond);
  riemann_solver_ = EnzoRiemann::construct_riemann(rsolver, cond);
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::setup_groups_(const EnzoFieldConditions cond)
{

  EnzoCenteredFieldRegistry registry;
  // all fields in primitive_group_ are permanent fields
  primitive_group_ = registry.build_cons_grouping(cond, "", "");
  // most (if not all) fields in conserved_group_ are temporary
  // For transparency, we create secondary "conserved" versions of fields that
  // represent quantities which are both conserved and primitive (e.g. density
  // and magnetic field). For example, primitive_group_ uses the permanent
  // field, "density", while conserved_group_ uses the temporary field
  // "cons_density"
  conserved_group_ = registry.build_cons_grouping(cond, "", "cons");

  bfieldi_group_ = new Grouping;
  bfieldi_group_->add("bfieldi_x", "bfield");
  bfieldi_group_->add("bfieldi_y", "bfield");
  bfieldi_group_->add("bfieldi_z", "bfield");

  // Get the names of the groups
  cons_group_names_ = registry.cons_group_names(cond, true);
  prim_group_names_ = registry.prim_group_names(cond, true);
}

//----------------------------------------------------------------------

EnzoMethodMHDVlct::~EnzoMethodMHDVlct()
{
  delete conserved_group_;
  delete bfieldi_group_;
  delete primitive_group_;
  delete half_dt_recon_;
  delete full_dt_recon_;
  delete riemann_solver_;
  delete eos_;
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  // a bug prevents us from pupping the groupings
  //p|conserved_group_;
  //p|bfieldi_group_;
  //p|primitive_group_;
  p|cond_;
  setup_groups_(cond_);

  p|eos_;
  p|half_dt_recon_;
  p|full_dt_recon_;
  p|riemann_solver_;

  p|cons_group_names_;
  p|prim_group_names_;
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute ( Block * block) throw()
{
  if (block->is_leaf()) {

    EnzoBlock * enzo_block = enzo::block(block);
    Field field = enzo_block->data()->field();

    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    // Make assertions about the size of the mesh
	   
    int ndim;
    if (enzo_block->GridDimension[2] == 1){
      ndim = 2;
    } else {
      ndim = 3;
    }
    

    // declaring Grouping that track temporary fields used for scratch space
    // Groupings of conserved fields and primitive fields are tracked as members
    // conserved_group_ and primitive_group_

    // left and right reconstructed primitives
    Grouping priml_group, primr_group;

    // flux fields
    Grouping xflux_group, yflux_group, zflux_group;

    // edge-centered electric fields
    Grouping efield_group;

    // Temporary field to store the central E-field (can reuse it to store
    // different components of the field at different times)
    std::string center_efield_name;

    // face-centered weight fields - Track the x/y/z upwind direction
    Grouping weight_group;

    // temp conserved group for storing values at the half time-step
    Grouping temp_conserved_group;

    // temp interface b-fileds
    Grouping temp_bfieldi_group;

    // left and right reconstructed conserved values
    // these are not strictly necessary
    Grouping consl_group, consr_group;

    // allocate the temporary fields (as necessary) and fill the field groupings
    allocate_temp_fields_(block, priml_group, primr_group, xflux_group,
			  yflux_group, zflux_group, efield_group,
			  center_efield_name, weight_group,
			  temp_conserved_group, temp_bfieldi_group,
			  consl_group, consr_group);

    // allocate constrained transport object
    EnzoConstrainedTransport ct = EnzoConstrainedTransport();

    EnzoFieldArrayFactory array_factory(block);

    // the following line is copied from EnzoMethodPpm & EnzoMethodHydro
    double dt = block->dt();

    eos_->conservative_from_primitive (block, *primitive_group_,
				       *conserved_group_);
    
    // Modify the following to work with groupings!
    // repeat the following twice (for half time-step and full time-step)
    for (int i=0;i<2;i++){
      double cur_dt;
      Grouping* cur_cons_group;
      Grouping* out_cons_group;
      Grouping* cur_bfieldi_group;
      Grouping* out_bfieldi_group;
      EnzoReconstructor *reconstructor;
      if (i == 0){
	cur_dt = dt/2.;
	out_cons_group = &temp_conserved_group;
	cur_cons_group = conserved_group_;
	out_bfieldi_group = &temp_bfieldi_group;
	cur_bfieldi_group = bfieldi_group_;
	reconstructor = half_dt_recon_;
      } else {
	cur_dt = dt;
	out_cons_group = conserved_group_;
	cur_cons_group = &temp_conserved_group;
	out_bfieldi_group = bfieldi_group_;
	cur_bfieldi_group = &temp_bfieldi_group;
	reconstructor = full_dt_recon_;
      }

      // Compute the primitive Quantites with Equation of State
      if (i == 1){
	eos_->primitive_from_conservative (block, *cur_cons_group,
					   *primitive_group_);
      }
      
      
      // Compute flux along each dimension
      compute_flux_(block, 0, *cur_cons_group, *cur_bfieldi_group, priml_group,
		    primr_group, xflux_group, consl_group, consr_group,
		    weight_group, *reconstructor);
      compute_flux_(block, 1, *cur_cons_group, *cur_bfieldi_group, priml_group,
		    primr_group, yflux_group, consl_group, consr_group,
		    weight_group, *reconstructor);
      if (ndim == 3){
	compute_flux_(block, 2, *cur_cons_group, *cur_bfieldi_group,
		      priml_group, primr_group, zflux_group, consl_group,
		      consr_group, weight_group, *reconstructor);
      }

      // Compute_efields
      //   - The first time, this implictly use the cell-centered primitives
      //     computed before the half timestep
      //   - The second time, it uses uses the the cell-centered primitives
      //     computed after the half timestep
      compute_efields_(block, xflux_group, yflux_group, zflux_group,
		       center_efield_name, efield_group, weight_group, ct);

      // Update quantities - add flux divergence
      update_quantities_(block, xflux_group, yflux_group, zflux_group,
			 *out_cons_group, cur_dt);

      // Add source terms - doesn't matter if done before/after adding flux
      // Update longitudinal B-field (add source terms of constrained transport)
      
      for (int dim = 0; dim<ndim; dim++){
	ct.update_bfield(block, dim, efield_group, *bfieldi_group_,
			 *out_bfieldi_group, cur_dt);
      }

      // Add other source terms?


      // Finally, update cell-centered B-field 
      for (int dim = 0; dim<ndim; dim++){
	ct.compute_center_bfield(block, dim, *out_cons_group,
				 *out_bfieldi_group);
      }
    }

    eos_->primitive_from_conservative (block, *conserved_group_,
				       *primitive_group_);

    //ASSERT("EnzoMethodMHDVlct","Early Exit",false);

    // Deallocate Temporary Fields
    deallocate_temp_fields_(block, priml_group, primr_group, xflux_group,
			    yflux_group, zflux_group, efield_group,
			    center_efield_name, weight_group,
			    temp_conserved_group, temp_bfieldi_group,
			    consl_group, consr_group);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::compute_flux_(Block *block, int dim,
				   Grouping &cur_cons_group,
				   Grouping &cur_bfieldi_group,
				   Grouping &priml_group,
				   Grouping &primr_group,
				   Grouping &flux_group,
				   Grouping &consl_group,
				   Grouping &consr_group,
				   Grouping &weight_group,
				   EnzoReconstructor &reconstructor)
{
  // First, reconstruct the left and right interface values
  reconstructor.reconstruct_interface(block, *primitive_group_, priml_group,
				      primr_group, dim, eos_);
  
  // Need to set the component of reconstructed B-field along dim, equal to
  // the corresponding longitudinal component of the B-field tracked at cell
  // interfaces (should probably be handled internally by reconstructor)
  EnzoFieldArrayFactory array_factory(block);
  EFlt3DArray temp,bfield, l_bfield,r_bfield;
  bfield = array_factory.interior_bfieldi(cur_bfieldi_group, dim);
  l_bfield = array_factory.reconstructed_field(priml_group, "bfield", dim, dim);
  r_bfield = array_factory.reconstructed_field(primr_group, "bfield", dim, dim);

  // All 3 array objects are the same shape
  // Iteration limits are generalized for 2D and 3D grids
  for (int iz=0; iz<bfield.shape(0); iz++) {
    for (int iy=0; iy<bfield.shape(1); iy++) {
      for (int ix=0; ix<bfield.shape(2); ix++) {
	l_bfield(iz,iy,ix) = bfield(iz,iy,ix);
	r_bfield(iz,iy,ix) = bfield(iz,iy,ix);
      }
    }
  }

  // Calculate reconstructed conserved values
  eos_->conservative_from_primitive(block,priml_group,consl_group);
  eos_->conservative_from_primitive(block,primr_group,consr_group);

  // Next, compute the fluxes
  riemann_solver_->solve(block, priml_group, primr_group, flux_group,
			 consl_group, consr_group, dim, eos_);
  
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

void EnzoMethodMHDVlct::compute_efields_(Block *block, Grouping &xflux_group,
				      Grouping &yflux_group,
				      Grouping &zflux_group,
				      std::string center_efield_name,
				      Grouping &efield_group,
				      Grouping &weight_group,
				      EnzoConstrainedTransport &ct)
{
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  int start_dim = 0;
  if (enzo_block->GridDimension[2] == 1) {
    // handling a 2D setup
    start_dim = 2;
  }

  // Maybe the following should be handled internally by ct?
  for (int i = start_dim; i < 3; i++){
    ct.compute_center_efield (block, i, center_efield_name, *primitive_group_);
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
			    weight_group);
  }
}

//----------------------------------------------------------------------

void EnzoMethodMHDVlct::update_quantities_(Block *block, Grouping &xflux_group,
					   Grouping &yflux_group,
					   Grouping &zflux_group,
					   Grouping &out_cons_group, double dt)
{
  // For now, not having density floor affect momentum or total energy density
  EnzoFieldArrayFactory array_factory(block);
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // cell-centered grid dimensions
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];

  const bool three_dim = (mz!=1);
  
  // Come up with iteration limits that are generalized for 2D grids
  int zstart, zstop;
  if (three_dim){
    zstart = 1;
    zstop = mz-1;
  } else {
    zstart = 0;
    zstop = mz;
  }

  // widths of cells
  enzo_float dtdx = dt/enzo_block->CellWidth[0];
  enzo_float dtdy = dt/enzo_block->CellWidth[1];
  enzo_float dtdz = dt/enzo_block->CellWidth[2];

  for (unsigned int group_ind=0; group_ind<cons_group_names_.size();
       group_ind++){
    // load group name and number of fields in the group
    std::string group_name = cons_group_names_[group_ind];
    if (group_name == "bfield"){
      continue;
    }
    int num_fields = conserved_group_->size(group_name);

    enzo_float floor_val =0;
    bool use_floor = false;
    if (group_name == "density"){
      floor_val = eos_->get_density_floor();
      use_floor=true;
    }

    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      // load in the quantities
      EFlt3DArray cur_cons, out_cons, xflux, yflux, zflux;
      cur_cons = array_factory.from_grouping(*conserved_group_, group_name,
					     field_ind);
      out_cons = array_factory.from_grouping(out_cons_group, group_name,
					     field_ind);
      xflux = array_factory.from_grouping(xflux_group, group_name, field_ind);
      yflux = array_factory.from_grouping(yflux_group, group_name, field_ind);
      if (three_dim){
	// only load if the grid is 3D
	zflux = array_factory.from_grouping(zflux_group, group_name, field_ind);
      }

      for (int iz=zstart; iz<zstop; iz++) {
	for (int iy=1; iy<my-1; iy++) {
	  for (int ix=1; ix<mx-1; ix++) {

	    enzo_float new_val;

	    if (three_dim){
	      new_val = (cur_cons(iz,iy,ix)
			 - dtdx * (xflux(iz,iy,ix) - xflux(iz,iy,ix-1))
			 - dtdy * (yflux(iz,iy,ix) - yflux(iz,iy-1,ix))
			 - dtdz * (zflux(iz,iy,ix) - zflux(iz-1,iy,ix)));
	    } else {
	      new_val = (cur_cons(iz,iy,ix)
			 - dtdx * (xflux(iz,iy,ix) - xflux(iz,iy,ix-1))
			 - dtdy * (yflux(iz,iy,ix) - yflux(iz,iy-1,ix)));
	    }
	    
	    if (use_floor){
	      new_val = std::max(new_val, floor_val);
	    }

	    out_cons(iz,iy,ix) = new_val;
	  }
	}
      }
    }
  }
  // Place floor on energy
  eos_->apply_floor_to_energy(block, out_cons_group);
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

// The temporary conserved fields have arbitrary initial values
// The calculation at the first half timestep fills in the inner cells in of
// the temporary fields. The outermost layer of cells on the grid are
// uninitialized. This can lead to NaNs or inf. To avoid this, the following
// helper function makes sure that the exterior layer has reasonable values.
// We could probably accomplish the same thing by setting the values to 0 or
// prim_floor. We also don't need to iterate over the full grid. None of this
// is strictly necessary since it will only be messing up the ghost zone - but
// if we want the calculation to be free of ever having nans/inf something like
// this is necessary.
//
// This is a quick solution - really only need to initalize values in the outer
// ghost zones
void copy_grouping_fields_(Block *block, Grouping &conserved_group,
			   Grouping &temp_cons_group,
			   std::vector<std::string> &group_names){

  EnzoFieldArrayFactory array_factory(block);
  
  for (unsigned int group_ind=0;group_ind<group_names.size();group_ind++){
    // load group name and number of fields in the group
    std::string group_name = group_names[group_ind];
    int num_fields = conserved_group.size(group_name);
     // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      // load in the quantities
      EFlt3DArray cur_cons, out_cons;
      cur_cons = array_factory.from_grouping(conserved_group, group_name,
					     field_ind);
      out_cons = array_factory.from_grouping(temp_cons_group, group_name,
					     field_ind);
      for (int iz=0; iz<cur_cons.shape(0); iz++) {
	for (int iy=0; iy<cur_cons.shape(1); iy++) {
	  for (int ix=0; ix<cur_cons.shape(2); ix++) {
	    out_cons(iz,iy,ix) = cur_cons(iz,iy,ix);
	  }
	}
      }
    }
  }
}

//----------------------------------------------------------------------

// helper function to initialize values of reconstructed primitives to floor
// that will never be initialized by reconstruction. This is necessary to avoid
// NaNs/inf because when EOS calculates reconstructed conserved quantites, EOS
// assumes that the reconstructed primitives are cell-centered (since they are
// not registered as face-centered) which means it will try to convert otherwise
// uninitialized values to conservatives.
//
// Really only need to initialize values at the end of the allocated blocks of
// memory.
void initialize_recon_prim_to_floor_(Block *block, Grouping &grouping,
				     EnzoEquationOfState &eos,
				     std::vector<std::string> &prim_group_names)
{
  int start = -1;
  int stop = -1;
  Field field = block->data()->field();

  for (unsigned int group_ind=0;group_ind<prim_group_names.size(); group_ind++){
    // load group name and number of fields in the group
    std::string group_name = prim_group_names[group_ind];
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

void EnzoMethodMHDVlct::allocate_temp_fields_(Block *block,
					      Grouping &priml_group,
					      Grouping &primr_group,
					      Grouping &xflux_group,
					      Grouping &yflux_group,
					      Grouping &zflux_group,
					      Grouping &efield_group,
					      std::string &center_efield_name,
					      Grouping &weight_group,
					      Grouping &temp_conserved_group,
					      Grouping &temp_bfieldi_group,
					      Grouping &consl_group,
					      Grouping &consr_group)
{
  std::vector<std::string> prim_group_names = prim_group_names_;
  
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();
  FieldDescr * field_descr = field.field_descr();
  
  // First, reserve/allocate temporary conserved-related fields
  // allocate applicable temporary fields in primitive_group_
  for (unsigned int i=0;i<cons_group_names_.size();i++){
    std::string group_name = cons_group_names_[i];
    int num_fields = conserved_group_->size(group_name);

    for (int j=0;j<num_fields;j++){
      // Determine field_name
      std::string field_name = conserved_group_->item(group_name,j);

      if (field.is_field(field_name) &&
	  field.is_permanent(field.field_id(field_name))){
	continue;
      }

      // allocate temporary field
      int ir_ = field_descr->insert_temporary(field_name);
      field.allocate_temporary(ir_);
    }
  }

  // Prepare the temporary conserved fields (used to store values at half dt)
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names_,
			    temp_conserved_group, "temp_",0,0,0);

  // Prepare temporary flux fields
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names_,
			    xflux_group, "xflux_",-1,0,0);
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names_,
			    yflux_group, "yflux_",0,-1,0);
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names_,
			    zflux_group, "zflux_",0,0,-1);

  // Prepare left and right reconstructed conserved fields
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names_,
			    consl_group, "cleft_",0,0,0);
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names_,
			    consr_group, "cright_",0,0,0);

  // Prepare temporary fields for priml and primr
  // As necessary, we pretend that these are centered along:
  //   - z and have shape (mz-1,  my,  mx)
  //   - y and have shape (  mz,my-1,  mx)
  //   - x and have shape (  mz,  my,mx-1)
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names_,
			    priml_group, "left_",0,0,0);
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names_,
			    primr_group, "right_",0,0,0);

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


  // Initialize Values (that would not otherwise be initialized) in order to
  // avoid NaN/Inf during calculations.

  // initialize values of the temporary conserved group (out in the ghost
  // zones to avoid NaNs during conversion to primitives  - these would
  // propogate through to the flux calculation)
  copy_grouping_fields_(block, *conserved_group_, temp_conserved_group,
			cons_group_names_);
  // initialize exterior face values of the temporary interface B-field
  // (this is just to avoid NaNs while computing the cell-centered B-fields
  // out in the ghost zone - this is NOT necessary if the underlying fields
  // are guarunteed to be set to 0 at start)
  std::vector<std::string> bfieldi_group_names{"bfield"};
  copy_grouping_fields_(block, *bfieldi_group_, temp_bfieldi_group,
			bfieldi_group_names);

  // initialize values of reconstructed primitives
  initialize_recon_prim_to_floor_(block, priml_group, *eos_, prim_group_names_);
  initialize_recon_prim_to_floor_(block, primr_group, *eos_, prim_group_names_);
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
						Grouping &xflux_group,
						Grouping &yflux_group,
						Grouping &zflux_group,
						Grouping &efield_group,
						std::string center_efield_name,
						Grouping &weight_group,
						Grouping &temp_conserved_group,
						Grouping &temp_bfieldi_group,
						Grouping &consl_group,
						Grouping &consr_group)
{

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // deallocate cell-centered conservative quantity fields
  deallocate_grouping_fields_(field, cons_group_names_, *conserved_group_);
  deallocate_grouping_fields_(field, cons_group_names_, temp_conserved_group);
  deallocate_grouping_fields_(field, cons_group_names_, xflux_group);
  deallocate_grouping_fields_(field, cons_group_names_, yflux_group);
  deallocate_grouping_fields_(field, cons_group_names_, zflux_group);
  deallocate_grouping_fields_(field, cons_group_names_, consl_group);
  deallocate_grouping_fields_(field, cons_group_names_, consr_group);

  // deallocate the (relevant) primitive primitive quantity fields
  deallocate_grouping_fields_(field, prim_group_names_, priml_group);
  deallocate_grouping_fields_(field, prim_group_names_, primr_group);

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
  // Implicitly assumes that "pressure" is a permanent field
  // analogous to ppm timestep calulation, probably want to require that cfast
  // is no smaller than some tiny positive number. 

  // Compute the pressure (assumes that "pressure" is a permanent field)
  //eos_->compute_pressure(block, *conserved_group_, *primitive_group_);
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

  pressure = array_factory.from_grouping(*primitive_group_, "pressure", 0);

  // Get iteration limits
  // Like ppm and ppml, access active region info from enzo_block attributes
  EnzoBlock * enzo_block = enzo::block(block);

  int ndim;
  if (enzo_block->GridDimension[2] == 1){
    ndim = 2;
  } else {
    ndim = 3;
  }

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

	// Using Gaussian units (where the magnetic permeability is mu=1)
	enzo_float inv_dens= 1./density(iz,iy,ix);
	enzo_float cfast = std::sqrt(gamma * pressure(iz,iy,ix) * inv_dens + 
				     bmag_sq * inv_dens);
 
	dtBaryons = std::min(dtBaryons,
			     dx/(std::fabs(velocity_x(iz,iy,ix))
					   + cfast));
	dtBaryons = std::min(dtBaryons,
			     dy/(std::fabs(velocity_y(iz,iy,ix))
					   + cfast));
	if (ndim == 3){
	  dtBaryons = std::min(dtBaryons,
			       dz/(std::fabs(velocity_z(iz,iy,ix))
					     + cfast));
	}
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
