
#include "cello.hpp"
#include "enzo.hpp"

#define TRACE_VLCT(MESSAGE)                                       \
 CkPrintf ("%s:%d TRACE_VLCT %s %s\n",                            \
           __FILE__,__LINE__, block->name().c_str(), MESSAGE);	  \
 fflush (stdout);

enzo_float* load_grouping_field_(Field *field, Grouping *grouping,
				 std::string group_name, int index)
{
  int size = grouping->size(group_name);
  if ((size == 0) || (size>=index)){
    return NULL;
  }
  return (enzo_float *) field->values(grouping->item(group_name,index));
}

// The idea is to completely replace the previous helper function with the
// following
void load_grouping_field_(Block *block, Grouping &grouping,
			  std::string group_name, int index,
			  EnzoArray<enzo_float> &array)
{
  Field field = block->data()->field();
  int size = grouping.size(group_name);
  if ((size == 0) || (size>=index)){
    return;
  }
  std::string field_name = grouping.item(group_name,index);
  enzo_float *data = (enzo_float *) field.values(field_name);
  int id = field.field_id(field_name);
  int mx, my, mz;
  field.dimensions (id,&mx,&my,&mz);
  array.initialize_wrapper(data,mz,my,mx);
}

// Because temporary fields must be allocated as cell-centered, we do not
// have FieldDescr internally track the centering of the fields. Consequently,
// we must specify the centering of the EnzoArray. Setting cell_centered_{dim}
// to true, it indicates that values are cell-centered along {dim}. Face
// -centered temporary fields do not have values at the exterior of the mesh.
void load_temp_interface_grouping_field_(Block *block, Grouping &grouping,
					 std::string group_name, int index,
					 EnzoArray<enzo_float> &array,
					 bool cell_centered_x,
					 bool cell_centered_y,
					 bool cell_centered_z)
{
  int size = grouping.size(group_name);
  if ((size == 0) || (size>=index)){
    return;
  }
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();
  enzo_float *data =(enzo_float *)field.values(grouping.item(group_name,index));
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];
  if (!cell_centered_x) {
    mx--;
  }
  if (!cell_centered_y) {
    my--;
  }
  if (!cell_centered_z) {
    mz--;
  }
  array.initialize_wrapper(data,mz,my,mx);
}

// Helper function designed to help initialize an EnzoArray representing a
// component of the interface (longitudinal) B-fields. 
void load_interior_bfieldi_field_(Block *block, Grouping &grouping,
				  int dim, EnzoArray<enzo_float> &array)
{
  int size = grouping.size("bfield");
  if ((size == 0) || (size>=dim)){
    return;
  }
  EnzoBlock * enzo_block = enzo::block(block);
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];
  int dix = (dim == 0) ? 1 : 0;
  int diy = (dim == 1) ? 1 : 0;
  int diz = (dix == diy) ? 1 : 0;

  Field field = enzo_block->data()->field();
  std::string field_name = grouping.item("bfield",dim);
  enzo_float *data =(enzo_float *)field.values(field_name);
  
  if (field.is_permanent(field.field_id(field_name))){
    EnzoArray<enzo_float> temp_array;
    temp_array.initialize_wrapper(data,mz+diz,my+diy,mx+dix);
    array.initialize_subarray(temp_array, diz, mz, diy, my, dix, mx);
  } else {
    array.initialize_wrapper(data, mz-diz, my-diy, mx-dix);
  }
}

// These 2 vectors need to be updated as more physics are added (e.g. dual
// energy formalism, CR Energy & Flux, species, colors)

// The following function sets up a list of cell-centered conserved groups
// conserved_group_, temp_conserved_group, xflux_group, yflux_group, and
// zflux_group will have the same subset of groups
std::vector<std::string> EnzoMethodVlct::cons_group_names={"density","momentum",
							   "total_energy",
							   "bfield"};

// The following function sets up a listof groups of primitive quantities
// primitive_group_, priml_group, and primr_group all will have the same
// subset of groups
std::vector<std::string> EnzoMethodVlct::prim_group_names={"density","velocity",
							   "pressure","bfield"};

//----------------------------------------------------------------------

EnzoMethodVlct::EnzoMethodVlct (double gamma)
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
  // total energy is the energy density
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
  bfieldi_group_ = new Grouping;
  primitive_group_ = new Grouping;
  // need a more sophisticated way to setup groups in the future
  setup_groups_();
  // Temporarilly use following default values
  double density_floor = 0;
  double pressure_floor = 0;
  
  eos_ = new EnzoEOSIdeal(gamma, density_floor, pressure_floor);
  half_dt_recon_ = new EnzoReconstructorPLM;
  full_dt_recon_ = new EnzoReconstructorPLM;
  riemann_solver_ = new EnzoRiemannHLLE;

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

  bfieldi_group_->add("bfield", "bfieldi_x");
  bfieldi_group_->add("bfield", "bfieldi_y");
  bfieldi_group_->add("bfield", "bfieldi_z");

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
  delete half_dt_recon_;
  delete full_dt_recon_;
  delete riemann_solver_;
  delete eos_;
}

//----------------------------------------------------------------------

void EnzoMethodVlct::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  // I think this is appropriate, but not totally sure
  // need to initialize to NULL
  p|eos_;
  p|half_dt_recon_;
  p|full_dt_recon_;
  p|riemann_solver_;

  // Currently rebuilding the groups every time. Need to pup conserved_group
  // and primitive_group_ in the future
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    EnzoBlock * enzo_block = enzo::block(block);
    int ndim;
    if (enzo_block->GridDimension[2] == 1){
      ndim = 2;
    } else {
      ndim = 3;
    }
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
			  center_efield_id, weight_group, temp_conserved_group,
			  temp_bfieldi_group, consl_group, consr_group);

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
      eos_->primitive_from_conservative (block, *cur_cons_group,
					 *primitive_group_);

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
		       center_efield_id, efield_group, weight_group, ct);

      // Update quantities - add flux divergence
      update_quantities_(block, xflux_group, yflux_group, zflux_group,
			 *out_cons_group, cur_dt);

      // Add source terms - doesn't matter if done before/after adding flux
      // Update longitudinal B-field (add source terms of constrained transport)
      for (int dim = 0; dim<ndim; dim++){
	ct.update_bfield(block, dim, efield_group, *bfieldi_group_,
			 *out_bfieldi_group, dt);
      }
      // Add other source terms?


      // Finally, update cell-centered B-field 
      for (int dim = 0; dim<ndim; dim++){
	ct.compute_center_bfield(block, dim, *out_cons_group,
				 *out_bfieldi_group, dt);
      }
    }

    // Deallocate Temporary Fields
    deallocate_temp_fields_(block, priml_group, primr_group, xflux_group,
			    yflux_group, zflux_group, efield_group,
			    center_efield_id, weight_group,
			    temp_conserved_group, temp_bfieldi_group,
			    consl_group, consr_group);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute_flux_(Block *block, int dim,
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

  EnzoArray<enzo_float> bfield, l_bfield,r_bfield;
  load_interior_bfieldi_field_(block, cur_bfieldi_group, dim, bfield);
  load_temp_interface_grouping_field_(block, priml_group, "bfield", dim,
				      l_bfield, dim == 0, dim == 1, dim == 2);
  load_temp_interface_grouping_field_(block, priml_group, "bfield", dim,
				      r_bfield, dim == 0, dim == 1, dim == 2);

  // All 3 array objects are the same shape
  // Iteration limits are generalized for 2D and 3D grids
  for (int iz=0; iz<bfield.length_dim2(); iz++) {
    for (int iy=0; iy<bfield.length_dim1(); iy++) {
      for (int ix=0; ix<bfield.length_dim0(); ix++) {
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

  EnzoArray<enzo_float> density_flux, weight_field;
  load_temp_interface_grouping_field_(block, flux_group, "density", 0,
				      density_flux, dim == 0, dim == 1,
				      dim == 2);
  load_temp_interface_grouping_field_(block, weight_group, "weight", dim,
				      weight_field, dim == 0, dim == 1,
				      dim == 2);

  // Iteration limits compatible with both 2D and 3D grids
  for (int iz=0; iz<density_flux.length_dim2(); iz++) {
    for (int iy=0; iy<density_flux.length_dim1(); iy++) {
      for (int ix=0; ix<density_flux.length_dim0(); ix++) {
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

void EnzoMethodVlct::compute_efields_(Block *block, Grouping &xflux_group,
				      Grouping &yflux_group,
				      Grouping &zflux_group,
				      int center_efield_id,
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
    ct.compute_center_efield (block, i, center_efield_id, *primitive_group_);
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

    ct.compute_edge_efield (block, i, center_efield_id, efield_group,
			    *jflux_group, *kflux_group, *primitive_group_,
			    weight_group);
  }
}

//----------------------------------------------------------------------

void EnzoMethodVlct::update_quantities_(Block *block, Grouping &xflux_group,
					Grouping &yflux_group,
					Grouping &zflux_group,
					Grouping &out_cons_group, double dt)
{
  // Need to address placing floors on density and energy_density
  std::vector<std::string> cons_group_names = EnzoMethodVlct::cons_group_names;

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // cell-centered iteration dimensions
  int mx = enzo_block->GridDimension[0];
  int my = enzo_block->GridDimension[1];
  int mz = enzo_block->GridDimension[2];

  // widths of cells
  enzo_float dtdx = dt/enzo_block->CellWidth[0];
  enzo_float dtdy = dt/enzo_block->CellWidth[1];
  enzo_float dtdz = dt/enzo_block->CellWidth[2];

  for (unsigned int group_ind=0;group_ind<cons_group_names.size();group_ind++){
    // load group name and number of fields in the group
    std::string group_name = prim_group_names[group_ind];
    if (group_name == "bfield"){
      continue;
    }
    int num_fields = conserved_group_->size(group_name);


    // iterate over the fields in the group
    for (int field_ind=0; field_ind<num_fields; field_ind++){

      // load in the quantities
      enzo_float *cur_cons = load_grouping_field_(&field, conserved_group_,
						  group_name, field_ind);
      enzo_float *out_cons = load_grouping_field_(&field, &out_cons_group,
						  group_name, field_ind);
      enzo_float *xflux = load_grouping_field_(&field, &xflux_group,
					       group_name, field_ind);
      enzo_float *yflux = load_grouping_field_(&field, &yflux_group,
					       group_name, field_ind);
      enzo_float *zflux = load_grouping_field_(&field, &zflux_group,
					       group_name, field_ind);

      for (int iz=1; iz<mz-1; iz++) {
	for (int iy=1; iy<my-1; iy++) {
	  for (int ix=1; ix<mx-1; ix++) {
	    int c_ind = ix + mx*(iy + my*iz);

	    out_cons[c_ind] = cur_cons[c_ind]
	      - dtdx * (xflux[ix+1 + (mx+1)*(iy + my*iz)]
			- xflux[ix + (mx+1)*(iy + my*iz)])
	      - dtdy * (yflux[ix + mx*(iy+1 + (my+1)*iz)]
			- yflux[ix + mx*(iy + (my+1)*iz)])
	      - dtdz * (zflux[ix + mx*(iy + my*(iz+1))]
			- zflux[ix + mx*(iy + my*iz)]);
	  }
	}
      }
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
			       std::vector<std::string> &group_names,
			       Grouping &grouping, std::string field_prefix)
{
  FieldDescr *field_descr = field.field_descr();
  for (unsigned int i=0;i<group_names.size();i++){
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
      field.allocate_temporary(ir_);

      // add the temporary field to grouping
      grouping.add(group_name,field_name);

    }
  }
}

// Helper function used to prepare groupings of temporary fields. These
// groupings only contain 3 fields - corresponding to vector components
//  - the names of the temporary fields are given by field_prefix + x,
//    field_prefix + y and field_prefix + z

void prep_temp_vector_grouping_(Field &field, std::string group_name,
				Grouping &grouping, std::string field_prefix)
{
  FieldDescr *field_descr = field.field_descr();
  for (int i=0;i<3;i++){
    // prepare field name
    std::string field_name;

    if (i == 0){
      field_name = field_prefix + "x";
    } else if (i == 1) {
      field_name = field_prefix + "y";
    } else {
      field_name = field_prefix + "z";
    }

    // Reserve a temporary field
    int ir_ = field_descr->insert_temporary(field_name);
    // Allocate the field
    field.allocate_temporary(ir_);

    grouping.add(group_name,field_name);
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
					   Grouping &weight_group,
					   Grouping &temp_conserved_group,
					   Grouping &temp_bfieldi_group,
					   Grouping &consl_group,
					   Grouping &consr_group)
{
  std::vector<std::string> cons_group_names = EnzoMethodVlct::cons_group_names;
  std::vector<std::string> prim_group_names = EnzoMethodVlct::prim_group_names;
  
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();
  FieldDescr * field_descr = field.field_descr();  

  // First, reserve/allocate temporary conserved-related fields

  // Prepare the temporary conserved fields (used to store values at half dt)
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    temp_conserved_group, "temp_");

  // Prepare temporary flux fields
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    xflux_group, "xflux_");
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    yflux_group, "yflux_");
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    zflux_group, "zflux_");

  // Prepare left and right reconstructed conserved fields
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    consl_group, "cleft_");
  prep_temp_field_grouping_(field, *conserved_group_, cons_group_names,
			    consr_group, "cright_");

  // Next, reserve/allocate temporary primitive-related fields
  // allocate applicable temporary fields in primitive_group_
  for (unsigned int i=0;i<prim_group_names.size();i++){
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
  // As necessary, we pretend that these are centered along:
  //   - z and have shape (mz-1,  my,  mx)
  //   - y and have shape (  mz,my-1,  mx)
  //   - x and have shape (  mz,  my,mx-1)
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names,
			    priml_group, "left_");
  prep_temp_field_grouping_(field, *primitive_group_, prim_group_names,
			    primr_group, "right_");

  // allocate temporary efield fields
  // reserve/allocate fields for edge-centered electric fields
  prep_temp_vector_grouping_(field, "efield", efield_group,
			     "temp_efield_");

  // reserve/allocate field for face-centered electric field
  center_efield_id = field_descr->insert_temporary();
  field.allocate_temporary(center_efield_id);


  // reserve/allocate fields for weight fields
  // these are face-centered fields that store the upwind/downwind direction
  prep_temp_vector_grouping_(field, "weight", weight_group, "temp_weight_");

  // reserve allocate temporary interface bfields fields
  prep_temp_vector_grouping_(field, "bfield", temp_bfieldi_group,
			     "temp_bfieldi_");
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
					     Grouping &weight_group,
					     Grouping &temp_conserved_group,
					     Grouping &temp_bfieldi_group,
					     Grouping &consl_group,
					     Grouping &consr_group)
{
  std::vector<std::string> cons_group_names= EnzoMethodVlct::cons_group_names;
  std::vector<std::string> prim_group_names= EnzoMethodVlct::prim_group_names;

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // deallocate cell-centered conservative quantity fields
  deallocate_grouping_fields_(field, cons_group_names, temp_conserved_group);
  deallocate_grouping_fields_(field, cons_group_names, xflux_group);
  deallocate_grouping_fields_(field, cons_group_names, yflux_group);
  deallocate_grouping_fields_(field, cons_group_names, zflux_group);
  deallocate_grouping_fields_(field, cons_group_names, consl_group);
  deallocate_grouping_fields_(field, cons_group_names, consr_group);

  // deallocate the (relevant) primitive primitive quantity fields
  deallocate_grouping_fields_(field, prim_group_names, *primitive_group_);
  deallocate_grouping_fields_(field, prim_group_names, priml_group);
  deallocate_grouping_fields_(field, prim_group_names, primr_group);

  // deallocate electric fields
  std::vector<std::string> efield_group_names{"efield"};
  deallocate_grouping_fields_(field, efield_group_names, efield_group);

  field.deallocate_temporary(center_efield_id);

  // deallocate weights
  std::vector<std::string> weight_group_names{"weight"};
  deallocate_grouping_fields_(field, weight_group_names, weight_group);

  // deallocate the temporary longitudinal bfields
  std::vector<std::string> bfieldi_group_names{"bfield"};
  deallocate_grouping_fields_(field, bfieldi_group_names, temp_bfieldi_group);
}

//----------------------------------------------------------------------

double EnzoMethodVlct::timestep ( Block * block ) const throw()
{
  TRACE_VLCT("timestep()");
  // Implicitly assumes that "pressure" is a permanent field
  // analogous to ppm timestep calulation, probably want to require that cfast
  // is no smaller than some tiny positive number. 

  // Compute the pressure (assumes that "pressure" is a permanent field)
  eos_->compute_pressure(block, *conserved_group_, *primitive_group_);
  enzo_float gamma = eos_->get_gamma();

  EnzoArray<enzo_float> density, momentum_x, momentum_y, momentum_z;
  EnzoArray<enzo_float> bfieldc_x, bfieldc_y, bfieldc_z;
  load_grouping_field_(block, *conserved_group_, "density", 0, density);
  load_grouping_field_(block, *conserved_group_, "momentum", 0, momentum_x);
  load_grouping_field_(block, *conserved_group_, "momentum", 1, momentum_y);
  load_grouping_field_(block, *conserved_group_, "momentum", 2, momentum_z);
  load_grouping_field_(block, *conserved_group_, "bfield", 0, bfieldc_x);
  load_grouping_field_(block, *conserved_group_, "bfield", 1, bfieldc_y);
  load_grouping_field_(block, *conserved_group_, "bfield", 2, bfieldc_z);

  EnzoArray<enzo_float> pressure;
  load_grouping_field_(block, *conserved_group_, "pressure", 0, pressure);

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

	// Using Gaussian units (where the magnetic permeability is mu=1)
	enzo_float inv_dens= 1./density(iz,iy,ix);
	enzo_float cfast = std::sqrt(gamma * pressure(iz,iy,ix) * inv_dens + 
				     bmag_sq * inv_dens);
 
	dtBaryons = std::min(dtBaryons,
			     dx/(std::fabs(inv_dens * momentum_x(iz,iy,ix)
					   + cfast)));
	dtBaryons = std::min(dtBaryons,
			     dy/(std::fabs(inv_dens * momentum_y(iz,iy,ix)
					   + cfast)));
	dtBaryons = std::min(dtBaryons,
			     dz/(std::fabs(inv_dens * momentum_z(iz,iy,ix)
					   + cfast)));
      }
    }
  }

  /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
  dtBaryons *= 0.5*courant_;
  //if (dtBaryon <= 0){
  //  ERROR2("EnzoMethodVlct's timestep method() yields an invalid value (%lf) "
  //	   "for %s\n", (double)dtBaryon, block->name().c_str());
  //}
  return dtBaryons;
}
