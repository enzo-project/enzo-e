
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

  // Temporarilly will set member variables to NULL
  eos_ = NULL;
  half_dt_recon_ = NULL;
  full_dt_recon_ = NULL;
  riemann_solver_ = NULL;
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
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    // declare vectors of the tracked conserved and face-centered B-fields:
    std::vector<int> cons_ids;
    std::vector<int> bface_ids;

    // Fill in cons_ids and bface_ids
    EnzoBlock * enzo_block = enzo::block(block);
    Field field = enzo_block->data()->field();

    // this might be an attribute worth storing:
    std::string cons_names[8] = {"density", "momentum_x", "momentum_y",
				 "momentum_z", "total_energy", "bfieldc_x",
				 "bfieldc_y", "bfieldc_z"};

    for (int i=0;i<8;i++){
	cons_ids.push_back( field.field_id(cons_names[i]) );
    }

    bface_ids.push_back(field.field_id("bfieldi_x"));
    bface_ids.push_back(field.field_id("bfieldi_y"));
    bface_ids.push_back(field.field_id("bfieldi_z"));


    // declaring vectors that track temporary field ids used for scratch space
    // Can be optimized by declaring size right off the bat
    // Not currently tracking arrays of field ids for 2 reasons:
    //   1. Transparency (More obvious which method uses which fields)
    //   2. Avoiding hassle of migrating attributes around
    
    // primitive cell-centered quantities
    std::vector<int> prim_ids;

    // left and right reconstructed primitives
    std::vector<int> priml_ids; std::vector<int> primr_ids;

    // flux ids
    std::vector<int> xflux_ids;
    std::vector<int> yflux_ids;
    std::vector<int> zflux_ids;

    // edge-centered electric fields
    std::vector<int> efield_ids;

    // Temporary field to store the central E-field (can reuse it to store
    // different components of the field at different times)
    int center_efield_id;

    // cell-centered conserved ids at the half time-step
    std::vector<int> temp_cons_ids;

    // face-centered B-fields at the half time-step
    std::vector<int> temp_bface_ids;

    // allocate the temporary fields (as necessary) and fill in the field_ids
    // need to update to allocate fluxes along each dimension!!!
    allocate_temp_fields_(block, prim_ids, priml_ids, primr_ids, xflux_ids,
			  yflux_ids, zflux_ids, efield_ids, center_efield_id,
			  temp_cons_ids, temp_bface_ids, cons_ids);

    // the following line is copied from EnzoMethodPpm & EnzoMethodHydro
    double dt = block->dt();

    // repeat the following twice (for half time-step and full time-step
    for (int i=0;i<2;i++){
      double cur_dt;
      std::vector<int> *out_cons_ids;
      std::vector<int> *out_bface_ids;
      if (i == 0){
	cur_dt = dt/2.;
	out_cons_ids = &temp_cons_ids;
	out_bface_ids = &temp_bface_ids;

	// By Default, prim_ids points to the density and cell-centered B-field
	// field ids stored in cons_ids - this is correct for first half ts

	// Compute the primitive Quantites with Equation of State
	eos_->primitive_from_conservative (block, cons_ids, prim_ids);

      } else {
	cur_dt = dt;
	out_cons_ids = &cons_ids;
	out_bface_ids = &bface_ids;

	// prim_ids needs to point to the density and cell-centered B-field
	// field ids stored in temp_cons_ids
	prim_ids[0] = temp_cons_ids[0]; // density
	for (int i=0;i<3;i++){
	  prim_ids[i+5] = temp_cons_ids[i+5]; // mag field
	}
	// Compute the primitive Quantites with Equation of State
	eos_->primitive_from_conservative (block, temp_cons_ids, prim_ids);
      }

      // Compute flux along each dimension
      compute_flux_(block, 0, prim_ids, bface_ids, priml_ids, primr_ids,
		    xflux_ids,half_dt_recon_);
      compute_flux_(block, 1, prim_ids, bface_ids, priml_ids, primr_ids,
		    yflux_ids,half_dt_recon_);
      compute_flux_(block, 2, prim_ids, bface_ids, priml_ids, primr_ids,
		    zflux_ids,half_dt_recon_);
    
      // Compute_efields
      compute_efields_(block, xflux_ids, yflux_ids, zflux_ids,
		       center_efield_id, efield_ids, prim_ids);
      // Possibly handle the CR source terms (and other hydro source terms)
      // right here
      
      // Update quantities
      update_quantities_(block, xflux_ids, yflux_ids, zflux_ids,
			 efield_ids, cons_ids, *out_cons_ids, bface_ids,
			 *out_bface_ids, cur_dt);
    }

    // Deallocate Temporary Fields
    deallocate_temp_fields_(block, prim_ids, priml_ids, primr_ids, xflux_ids,
			    yflux_ids, zflux_ids, efield_ids, center_efield_id,
			    temp_cons_ids, temp_bface_ids);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute_flux_(Block *block, int dim,
				   const std::vector<int> &prim_ids,
				   const std::vector<int> &bface_ids,
				   std::vector<int> &priml_ids,
				   std::vector<int> &primr_ids,
				   std::vector<int> &flux_ids,
				   EnzoReconstructor *reconstructor)
{
  // To be filled in

  // we need to construct a list 
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute_efields_(Block *block,
				      const std::vector<int> &xflux_ids,
				      const std::vector<int> &yflux_ids,
				      const std::vector<int> &zflux_ids,
				      int center_efield_id,
				      const std::vector<int> &efield_ids,
				      const std::vector<int> &prim_ids)
{
  // To be filled in
}

//----------------------------------------------------------------------

void EnzoMethodVlct::update_quantities_(Block *block,
					const std::vector<int> &xflux_ids,
					const std::vector<int> &yflux_ids,
					const std::vector<int> &zflux_ids,
					const std::vector<int> &efield_ids,
					const std::vector<int> &cur_cons_ids,
					const std::vector<int> &out_cons_ids,
					const std::vector<int> &cur_bface_ids,
					const std::vector<int> &out_bface_ids,
					double dt)
{
  // To be filled in
  // Should probably pass in the reference conserved ids
  // Also pass in the output conserved ids (can optionally be the same as above)
}

//----------------------------------------------------------------------

void EnzoMethodVlct::allocate_temp_fields_(Block *block,
					   std::vector<int> &prim_ids,
					   std::vector<int> &priml_ids,
					   std::vector<int> &primr_ids,
					   std::vector<int> &xflux_ids,
					   std::vector<int> &yflux_ids,
					   std::vector<int> &zflux_ids,
					   std::vector<int> &efield_ids,
					   int &center_efield_id,
					   std::vector<int> &temp_cons_ids,
					   std::vector<int> &temp_bface_ids,
					   const std::vector<int> &cons_ids)
{
  // This will need to be modified when we add on more physics (e.g. CRs)

  FieldDescr * field_descr = cello::field_descr();
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  // First, reserve/allocate temporary conserved cell-centered fields
  for (int i=0;i<8;i++){
    // Reserve a temporary field
    int ir_ = field_descr->insert_temporary();
    // Allocate the field
    field.allocate_temporary(ir_);
    temp_cons_ids.push_back(ir_);
  }


  // Next, reserve/allocate temporary Face-Centered B-fields fields
  for (int i=0;i<3;i++){
    // Reserve a temporary field
    int ir_ = field_descr->insert_temporary();
    // Specify that the field is face-centered
    int cx = (i == 0) ? 1 : 0;
    int cy = (i == 1) ? 1 : 0;
    int cz = (i == 2) ? 1 : 0;
    field_descr->set_centering(ir_,cx,cy,cz);
    // Allocate the field
    field.allocate_temporary(ir_);
    temp_bface_ids.push_back(ir_);
  }

  // deal with cell-centered primitives
  // reuse density and cell-centered B-fields from cons_ids
  // need to take special care to swap density and cell-centered B-field
  // field ids after the half time step.
  for (int i=0;i<8;i++){
    if ((i > 0) && (i<5)) {
      // Reserve a temporary field
      int ir_ = field_descr->insert_temporary();
      // Allocate the field
      field.allocate_temporary(ir_);
      prim_ids.push_back(ir_);
    } else {
      // set the primitive id equal to that of the corresponding field
      prim_ids.push_back(cons_ids[i]);
    }
  }

  // reserve/allocate fields for reconstructed left/right interface primitives
  // we "cheat" here and say fields are corner-centered (so we can reuse the
  // temporary fields for each dimension)

  for (int i=0;i<2;i++){
    for (int j=0;i<8;j++){
      // Reserve a temporary field
      int ir_ = field_descr->insert_temporary();
      // Specify that the field is corner-centered
      field_descr->set_centering(ir_,1,1,1);
      // Allocate the field
      field.allocate_temporary(ir_);
      switch(i) {
      case 0 : priml_ids.push_back(ir_);
      case 1 : primr_ids.push_back(ir_);
      }
    }
  }

  // reserve/allocate fields for reconstructed for fluxes
  for (int i=0;i<3;i++){
    int cx = (i == 0) ? 1 : 0;
    int cy = (i == 1) ? 1 : 0;
    int cz = (i == 2) ? 1 : 0;
    for (int j=0;j<8;j++){
      // Reserve a temporary field
      int ir_ = field_descr->insert_temporary();
      // Specify that the field is face-centered
      field_descr->set_centering(ir_,cx,cy,cz);
      // Allocate the field
      field.allocate_temporary(ir_);

      switch(i) {
      case 0 : xflux_ids.push_back(ir_);
      case 1 : yflux_ids.push_back(ir_);
      case 2 : zflux_ids.push_back(ir_);
      }
    }
  }

  // reserve/allocate fields for edge-centered electric fields
  for (int i=0;i<3;i++){
    // Reserve a temporary field
    int ir_ = field_descr->insert_temporary();
    // Specify that the field is edge-centered
    int cx = (i == 0) ? 0 : 1;
    int cy = (i == 1) ? 0 : 1;
    int cz = (i == 2) ? 0 : 1;
    field_descr->set_centering(ir_,cx,cy,cz);
    // Allocate the field
    field.allocate_temporary(ir_);
    temp_bface_ids.push_back(ir_);
  }

  // reserve/allocate field for face-centered electric field
  center_efield_id = field_descr->insert_temporary();
  field.allocate_temporary(center_efield_id);
}

//----------------------------------------------------------------------

void EnzoMethodVlct::deallocate_temp_fields_(Block *block,
					     std::vector<int> &prim_ids,
					     std::vector<int> &priml_ids,
					     std::vector<int> &primr_ids,
					     std::vector<int> &xflux_ids,
					     std::vector<int> &yflux_ids,
					     std::vector<int> &zflux_ids,
					     std::vector<int> &efield_ids,
					     int &center_efield_id,
					     std::vector<int> &temp_cons_ids,
					     std::vector<int> &temp_bface_ids)
{
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  for (int i=0;i<8;i++){
    field.deallocate_temporary(temp_cons_ids[i]);
  }
  temp_cons_ids.clear();

  for (int i=0;i<3;i++){
    field.deallocate_temporary(temp_bface_ids[i]);
  }
  temp_bface_ids.clear();

  
  // need to be careful with prim_ids
  // we already deallocated the other temporary fields with temp_cons_ids
  for (int i=1;i<5;i++){
    field.deallocate_temporary(prim_ids[i]);
  }
  prim_ids.clear();

  for (int i=0;i<8;i++){
    field.deallocate_temporary(priml_ids[i]);
  }
  priml_ids.clear();

  for (int i=0;i<8;i++){
    field.deallocate_temporary(primr_ids[i]);
  }
  primr_ids.clear();

  for (int i=0;i<3;i++){
    std::vector<int> *flux_ids;
    switch(i){
    case 0 : flux_ids = &xflux_ids;
    case 1 : flux_ids = &yflux_ids;
    case 2 : flux_ids = &zflux_ids;
    }
    for (int i=0;i<8;i++){
      field.deallocate_temporary((*flux_ids)[i]);
    }
    flux_ids->clear();
  }

  for (int i=0;i<3;i++){
    field.deallocate_temporary(efield_ids[i]);
  }
  efield_ids.clear();

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
