
#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodVlct::EnzoMethodVlct ()
  : Method()
{
  // Initialize the default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
  			       enzo_sync_id_method_vlct);

  // Not clear if all fields should be refreshed (like in PPML constructor)
  FieldDescr * field_descr = cello::field_descr();
  
  refresh(ir)->add_field(field_descr->field_id("density"));
  refresh(ir)->add_field(field_descr->field_id("velocity_x"));
  refresh(ir)->add_field(field_descr->field_id("velocity_y"));
  refresh(ir)->add_field(field_descr->field_id("velocity_z"));
  // total energy should be specific energy
  refresh(ir)->add_field(field_descr->field_id("total_energy"));
  refresh(ir)->add_field(field_descr->field_id("bfieldc_x"));
  refresh(ir)->add_field(field_descr->field_id("bfieldc_y"));
  refresh(ir)->add_field(field_descr->field_id("bfieldc_z"));
  refresh(ir)->add_field(field_descr->field_id("bfieldi_x"));
  refresh(ir)->add_field(field_descr->field_id("bfieldi_y"));
  refresh(ir)->add_field(field_descr->field_id("bfieldi_z"));
  refresh(ir)->add_field(field_descr->field_id("acceleration_x"));
  refresh(ir)->add_field(field_descr->field_id("acceleration_y"));
  refresh(ir)->add_field(field_descr->field_id("acceleration_z"));

  conserved_group_ = new Grouping;
  primitive_group_ = new Grouping;
  // need a more sophisticated way to setup groups in the future
  setup_groups_();
  // Temporarilly will set member variables to NULL
  eos_ = NULL;
  half_dt_recon_ = NULL;
  full_dt_recon_ = NULL;
  riemann_solver_ = NULL;
}

//----------------------------------------------------------------------

void EnzoMethodVlct::setup_groups_()
{
  // Fill in entries in conserved_group_ and primitive_group_
  // conserved_group_ will be entirely permanent fields while primitive_group_
  // will be partially made of temporary fields

  // In the future, this may also incorporate species, colors, and CRs

  conserved_group_->add("density", "density");
  conserved_group_->add("momentum", "momentum_x");
  conserved_group_->add("momentum", "momentum_y");
  conserved_group_->add("momentum", "momentum_z");
  conserved_group_->add("total_energy", "total_energy");
  conserved_group_->add("bfield", "bfieldc_x");
  conserved_group_->add("bfield", "bfieldc_y");
  conserved_group_->add("bfield", "bfieldc_y");

  conserved_group_->add("bfieldi", "bfieldi_x");
  conserved_group_->add("bfieldi", "bfieldi_y");
  conserved_group_->add("bfieldi", "bfieldi_y");

  // For transparency we will create secondary "primitive" density and magnetic
  // fields
  primitive_group_->add("density", "prim_density");
  primitive_group_->add("velocity", "velocity_x");
  primitive_group_->add("velocity", "velocity_y");
  primitive_group_->add("velocity", "velocity_z");
  primitive_group_->add("pressure", "pressure");
  // The following are cell-centered
  primitive_group_->add("bfield", "prim_bfield_x");
  primitive_group_->add("bfield", "prim_bfield_y");
  primitive_group_->add("bfield", "prim_bfield_y");

  
}

//----------------------------------------------------------------------

EnzoMethodVlct::~EnzoMethodVlct()
{
  delete conserved_group_;
  delete primitive_group_;
}

//----------------------------------------------------------------------

void EnzoMethodVlct::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  // I think this is appropriate, but not totally sure
  // need to initialize to NULL
  //p|eos_;
  //p|half_dt_recon_;
  //p|full_dt_recon_;
  //p|riemann_solver_;

  // Currently rebuilding the groups every time. Need to pup conserved_group
  // and primitive_group_ in the future
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    EnzoBlock * enzo_block = enzo::block(block);
    Field field = enzo_block->data()->field();

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
    int center_efield_id;

    // temp conserved group for storing values at the half time-step
    Grouping* temp_conserved_group;

    // allocate the temporary fields (as necessary) and fill the field groupings
    allocate_temp_fields_(block, priml_group, primr_group, xflux_group,
			  yflux_group, zflux_group, efield_group,
			  center_efield_id, temp_conserved_group);

    // allocate constrained transport object
    EnzoConstrainedTransport ct = EnzoConstrainedTransport();
    
    // the following line is copied from EnzoMethodPpm & EnzoMethodHydro
    double dt = block->dt();

    // Modify the following to work with groupings!
    // repeat the following twice (for half time-step and full time-step)
    for (int i=0;i<2;i++){
      double cur_dt;
      Grouping* cur_cons_group;
      Grouping* out_cons_group;
      EnzoReconstructor *reconstructor;
      if (i == 0){
	cur_dt = dt/2.;
	out_cons_group = &temp_conserved_group;
	cur_cons_group = conserved_group_;
	reconstructor = half_dt_recon_;
      } else {
	cur_dt = dt;
	out_cons_group = conserved_group_;
	cur_cons_group = temp_conserved_group;
	reconstructor = full_dt_recon_;
      }

      // Compute the primitive Quantites with Equation of State
      eos_->primitive_from_conservative (block, *cur_cons_group,
					 *primitive_group_);

      // Compute flux along each dimension
      compute_flux_(block, 0, *cur_cons_group, priml_group, primr_group,
		    xflux_group, *reconstructor);
      compute_flux_(block, 1, *cur_cons_group, priml_group, primr_group,
		    yflux_group, *reconstructor);
      compute_flux_(block, 2, *cur_cons_group, priml_group, primr_group,
		    zflux_group, *reconstructor);

      // Compute_efields
      compute_efields_(block, xflux_group, yflux_group, zflux_group,
		       center_efield_id, efield_group, ct);

      // Update quantities
      update_quantities_(block, xflux_group, yflux_group, zflux_group,
			 *out_cons_group, cur_dt);

      // Add source Terms
    }

    // Deallocate Temporary Fields
    deallocate_temp_fields_(block, priml_group, primr_group, xflux_group,
			    yflux_group, zflux_group, efield_group,
			    center_efield_id, temp_conserved_group);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute_flux_(Block *block, int dim,
				   Grouping &cur_cons_group,
				   Grouping &priml_group,
				   Grouping &primr_group,
				   Grouping &flux_group,
				   EnzoReconstructor &reconstructor)
{
  // First, let's reconstruct the left and right interface values
  reconstructor.reconstruct_interface(block, *primitive_group_, priml_group,
				      primr_group, dim);

  // Need to set the component of reconstructed B-field along dim, equal to
  // the corresponding longitudinal component of the B-field tracked at cell
  // interfaces
  // This should probably be handled internally by reconstructor
  // NEED TO CHECK THAT THE DIMENSIONS GIVEN FOR THE FLUX IDS ARE TRULY
  // FACE-CENTERED

  Field field = block->data()->field();
  const int id = field.field_id(cur_cons_group.item("bfieldi",dim));;

  // iteration dimensions (includes ghost zones)
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);

  // Get the field data
  enzo_float *bfield = load_grouping_field_(&field, &cur_cons_group, "bfieldi",
					    dim);
  enzo_float *l_bfield = load_grouping_field_(&field, &priml_group, "bfield",
					      dim);
  enzo_float *r_bfield = load_grouping_field_(&field, &primr_group, "bfield",
					      dim);

  for (int iz=1; iz<mz-1; iz++) {
    for (int iy=1; iy<my-1; iy++) {
      for (int ix=1; ix<mx-1; ix++) {
	int i = ix + mx*(iy + my*iz);
	l_bfield[i] = bfield[i];
	r_bfield[i] = bfield[i];
      }
    }
  }

  // Next, compute the fluxes
  riemann_solver_->solve(block, priml_group, primr_group, flux_group, dim,
			 eos_);
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute_efields_(Block *block, Grouping &xflux_group,
				      Grouping &yflux_group,
				      Grouping &zflux_group,
				      int center_efield_id,
				      Grouping &efield_group,
				      EnzoConstrainedTransport &ct)
{
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // Maybe the following should be handled internally by ct?
  for (int i = 0; i <3; i++){
    ct.compute_cell_center_efield (block, i, center_efield_id,
				   *primitive_group_);
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

    int efield_id = field.field_id(efield_group.item("efield",i));

    ct.compute_edge_efield (block, i, efield_id, center_efield_id, *jflux_group,
			    *kflux_group, *primitive_group_);
  }
}

//----------------------------------------------------------------------

void EnzoMethodVlct::update_quantities_(Block *block, Grouping &xflux_group,
					Grouping &yflux_group,
					Grouping &zflux_group,
					Grouping &out_cons_group, double dt)
{
  // To be filled in
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
			       std::string *group_names, int num_groups,
			       Grouping &grouping, std::string field_prefix,
			       int cx, int cy, int cz)
{
  for (int i=0;i<num_groups;i++){
    std::string group_name = group_names[i];
    int num_fields = ref_grouping.size(group_name);
    if (num_fields == 0) {
      continue; // This group is not tracked by ref_grouping
    }

    for (int j=0;j<num_fields;j++){
      // Determine field_name
      std::string field_name = field_prefix + (ref_grouping.item(group_name,j));

      // allocate temporary field
      int ir_ = field_descr->insert_temporary(field_name);
      field_descr->set_centering(ir_,cx,cy,cz);
      field.allocate_temporary(ir_);

      // add the temporary field to grouping
      grouping[j].add(group,field_name);

    }
  }
}

void EnzoMethodVlct::allocate_temp_fields_(Block *block,
					   Grouping &priml_group,
					   Grouping &primr_group,
					   Grouping &xflux_group,
					   Grouping &yflux_group,
					   Grouping &zflux_group,
					   Grouping &efield_group,
					   int &center_efield_id,
					   Grouping &temp_conserved_group);
{
  // This will need to be modified when we add on more physics (e.g. CRs)

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();
  FieldDescr * field_descr = field.field_descr();  

  // First, reserve/allocate temporary conserved-related fields

  // Prepare the temporary conserved fields (used to store values at half dt)
  // This is the only grouping other that conserved_group_ to include both
  // face-centered and cell-centered bfields
  // --Prepare cell-centered fields:
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    num_prim_group_names, temp_conserved_group,
			    "temp_", 0, 0, 0);

  // --Prepare interface B-fields:
  for (int i=0;i<3;i++){
    // Determine Field name
    std::string field_name = "temp_" + (conserved_group_->item("bfieldi",i));

    // Reserve a temporary field
    int ir_ = field_descr->insert_temporary(field_name);

    // Specify that the field is face-centered
    int cx = (i == 0) ? 1 : 0;
    int cy = (i == 1) ? 1 : 0;
    int cz = (i == 2) ? 1 : 0;
    field_descr->set_centering(ir_,cx,cy,cz);

    // Allocate the field
    field.allocate_temporary(ir_);

    // Add field to group
    temp_conserved_group.add("bfieldi",field_name);
  }


  // Prepare temporary flux fields [They should not include interface bfields]
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    num_prim_group_names, xflux_group, "xflux_",
			    1, 0, 0);
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    num_prim_group_names, yflux_group, "yflux_",
			    0, 1, 0);
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    num_prim_group_names, zflux_group, "zflux_",
			    0, 0, 1);
  
  // allocate applicable temporary fields in primitive_group_
  for (int i=0;i<num_prim_group_names;i++){
    std::string group_name = prim_group_names[i];
    int num_fields = primitive_group_->size(group_name);
    if (num_fields == 0) {
      continue; // This group is not tracked by ref_grouping
    }

    for (int j=0;j<num_fields;j++){
      // Determine field_name
      std::string field_name = primitive_group_->item(group_name,j);

      // slightly unsure if the field_id will return -1 if the temporary field
      // was used and deallocated in a prior timestep
      if (field.field_id(field_name) != -1){
	continue;
      }

      // allocate temporary field
      int ir_ = field_descr->insert_temporary(field_name);
      field.allocate_temporary(ir_);
    }
  }


  // Prepare temporary fields for priml and primr
  // We "cheat" here and make them corner-centered so we can be sure that
  // there is memory to reuse the fields and treat them as face-centered in all
  // directions
  prep_temp_field_grouping_(field, *primitive_group_, temp_prim_group_names,
			    num_prim_group_names, priml_group, "left_",
			    1, 1, 1);
  prep_temp_field_grouping_(field, *primitive_group_, temp_prim_group_names,
			    num_prim_group_names, primr_group, "right_",
			    1, 1, 1);

  // allocate temporary efield fields
  // Fluxes, reconstructed left/right interface primitives, fluxes, and efields
  // can be better optimized. Currently, they are constructed so
  // that they will have the exact memory layout expected for a face-centered
  // field (including the ghost zones). We are using this for right now, for
  // simplicity. (This is primarily useful for going between face-centered
  // B fields and these temporary fields).
  //
  // In fact this is not necessary. When we are computing values at the
  // interfaces between cells, we only need N-1 elements along that dimension
  // (we don't calculate anything for outermost faces)

  // reserve/allocate fields for edge-centered electric fields
  for (int i=0;i<3;i++){
    // get field_name and prepare edge-centered info
    std::string field_name;
    int cx=1; int cy=1; int cz = 1;
    if (i == 0){
      cx = 0;
      field_name = "temp_efield_x";
    } else if (i == 1) {
      cy = 0;
      field_name = "temp_efield_y";
    } else {
      cz = 0;
      field_name = "temp_efield_z";
    }

    // Reserve a temporary field
    int ir_ = field_descr->insert_temporary(field_name);
    field_descr->set_centering(ir_,cx,cy,cz);
    // Allocate the field
    field.allocate_temporary(ir_);

    efield_group.add("efieldi",field_name);
  }

  // reserve/allocate field for face-centered electric field
  center_efield_id = field_descr->insert_temporary();
  field.allocate_temporary(center_efield_id);
}

//----------------------------------------------------------------------

// Helper function that deallocates the temporary fields of a grouping
void deallocate_grouping_fields_(Field &field, std::string *group_names,
				 int num_groups, Grouping &grouping)
{
  for (int i=0;i<num_groups;i++){
    std::string group_name = group_names[i];
    int num_fields = grouping.size(group_name);
    if (num_fields == 0) {
      continue; // This group is not tracked by grouping
    }

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

void EnzoMethodVlct::deallocate_temp_fields_(Block *block,
					     Grouping &priml_group,
					     Grouping &primr_group,
					     Grouping &xflux_group,
					     Grouping &yflux_group,
					     Grouping &zflux_group,
					     Grouping &efield_group,
					     int center_efield_id,
					     Grouping &temp_conserved_group);
{
  // Need to update
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // deallocate cell-centered conservative quantity fields
  deallocate_grouping_fields_(field, cons_group_names, num_cons_groups,
			      temp_conserved_group);
  deallocate_grouping_fields_(field, cons_group_names, num_cons_groups,
			      xflux_group);
  deallocate_grouping_fields_(field, cons_group_names, num_cons_groups,
			      yflux_group);
  deallocate_grouping_fields_(field, cons_group_names, num_cons_groups,
			      zflux_group);

  // deallocate the temporary longitudinal bfields
  std::string bfield_group_names[1] = {"bfieldi"};
  deallocate_grouping_fields_(field, bfield_group_names, 1,
			      temp_conserved_group);

  // deallocate the (relevant) primitive primitive quantity fields
  deallocate_grouping_fields_(field, prim_group_names, num_prim_groups,
			      primitive_group_);
  deallocate_grouping_fields_(field, prim_group_names, num_prim_groups,
			      priml_group);
  deallocate_grouping_fields_(field, prim_group_names, num_prim_groups,
			      primr_group);

  // deallocate electric fields
  std::string efield_group_names[1] = {"efield"};
  deallocate_grouping_fields_(field, efield_group_names, 1, efield_group);

  field.deallocate_temporary(center_efield_id);
}

//----------------------------------------------------------------------

double EnzoMethodVlct::timestep ( Block * block ) const throw()
{
  /* Points to address:
   *  1. Make sure the overloading of std::sqrt works correctly [and always 
   *     returns the precision of enzo_float
   *  2. Possibly rename the fields 
   *  3. Guarantee that cfast is positive */

  // initialize

  enzo_float dtBaryons = ENZO_HUGE_VAL;

  /* Compute the pressure. */

  EnzoBlock * enzo_block = enzo::block(block);

  const int in = cello::index_static();

  enzo_float gamma = EnzoBlock::Gamma[in];
  EnzoComputePressure compute_pressure (gamma, false);
  compute_pressure.compute(enzo_block);

  Field field = enzo_block->data()->field();

  enzo_float * density    = (enzo_float *) field.values("density");
  enzo_float * velocity_x = (enzo_float *) field.values("velocity_x");
  enzo_float * velocity_y = (enzo_float *) field.values("velocity_y");
  enzo_float * velocity_z = (enzo_float *) field.values("velocity_z");
  enzo_float * bfieldc_x = (enzo_float *) field.values("bfieldc_x");
  enzo_float * bfieldc_y = (enzo_float *) field.values("bfieldc_y");
  enzo_float * bfieldc_z = (enzo_float *) field.values("bfieldc_z");
  enzo_float * pressure = (enzo_float *) field.values("pressure");

  /* Like ppm and ppml, access active region info from enzo_block attributes
   * In these attributes 0,1,2 correspond to x,y,z */

  /* x and y dimensions */
  int xdim, ydim;
  xdim = enzo_block->GridDimension[0]; ydim = enzo_block->GridDimension[1];

  /* The start index of the active region */
  int ix_start, iy_start, iz_start;
  ix_start = enzo_block->GridStartIndex[0];
  iy_start = enzo_block->GridStartIndex[1];
  iz_start = enzo_block->GridStartIndex[2];

  /* The end index (final valid index) of the active region */
  int ix_end, iy_end, iz_end;
  ix_end = enzo_block->GridEndIndex[0];
  iy_end = enzo_block->GridEndIndex[1];
  iz_end = enzo_block->GridEndIndex[2];

  /* widths of cells */
  double dx,dy,dz;
  dx = enzo_block->CellWidth[0];
  dy = enzo_block->CellWidth[1];
  dz = enzo_block->CellWidth[2];

  /* timestep is the minimum of 0.5 * dw/(abs(vw)+cfast) for all dimensions.
   * hw and vw are the the width of the cell and velocity along dimension w.
   * cfast = fast magnetosonic speed [Convention is to use max value: 
   * cfast = (va^2+cs^2) ] */

  for (int iz=iz_start; iz<=iz_end; iz++) {
    for (int iy=iy_start; iy<=iy_end; iy++) {
      for (int ix=ix_start; ix<=ix_end; ix++) {
	int i = ix + xdim*(iy + ydim*iz);
	enzo_float bmag_sq = (bfieldc_x[i] * bfieldc_x[i] +
			      bfieldc_y[i] * bfieldc_y[i] +
			      bfieldc_z[i] * bfieldc_z[i]);

	/* analogous to ppm timestep calulation, probably want to require that
	 * cfast is no smaller than some tiny positive number. */
	enzo_float cfast = std::sqrt(gamma * pressure[i] / density[i]  + 
				     bmag_sq /(4 * cello::pi * density[i]));

	dtBaryons = std::min(dtBaryons, dx/(std::fabs(velocity_x[i]) + cfast));
	dtBaryons = std::min(dtBaryons, dy/(std::fabs(velocity_y[i]) + cfast));
	dtBaryons = std::min(dtBaryons, dz/(std::fabs(velocity_z[i]) + cfast));
      }
    }
  }

  /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
  dtBaryons *= 0.5*courant_;
  return dtBaryons;
}
