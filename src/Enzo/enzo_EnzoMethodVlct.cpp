
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

    
    // To be filled in!

    

    

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
    std::vector<int> flux_ids;

    // edge-centered electric fields
    std::vector<int> efield_ids;

    // Temporary field to store the central E-field (can reuse it to store
    // different components of the field at different times)
    int center_efield_id;
    // conserved ids at the half time-step
    std::vector<int> temp_cons_ids;

    // allocate the temporary fields (as necessary) and fill in the field_ids
    allocate_temp_fields_(block, prim_ids, priml_ids, primr_ids, flux_ids,
			  efield_ids, center_efield_id, temp_cons_ids);
    
    // the tracked conserved ids:
    std::vector<int> cons_ids;

    // Fill in the conserved ids

    // Compute the primitive quantities with the EOS


    // First Half-ts:
    //   Compute the Conserved Quantites with Equation of State
    //   Then, for each dimension, call compute_flux_ method, for each dimension
    //   to compute the fluxes
    //   After that, compute_efields
    //   Update quantities
    //   Possibly handle the CRs

    // Repeat the operations of the First Half-ts, but now for the full timestep

    // Deallocate Temporary Fields
    deallocate_temp_fields_(block, prim_ids, priml_ids, primr_ids, flux_ids,
			    efield_ids, center_efield_id, temp_cons_ids);
  }

  block->compute_done();
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute_flux_(Block *block, int flux_id, int dim,
				   const std::vector<int> &prim_ids,
				   EnzoReconstructor *reconstructor)
{
  // To be filled in

  // we need to construct a list 
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute_efields_(Block *block,
				      const std::vector<int> &flux_ids,
				      int center_efield_id,
				      const std::vector<int> &efield_ids,
				      const std::vector<int> &prim_ids)
{
  // To be filled in
}

//----------------------------------------------------------------------

void EnzoMethodVlct::update_quantities_(Block *block,
					const std::vector<int> &flux_ids,
					const std::vector<int> &efield_ids,
					const std::vector<int> &cur_cons_ids,
					const std::vector<int> &out_cons_ids,
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
					   std::vector<int> &flux_ids,
					   std::vector<int> &efield_ids,
					   int &center_efield_id,
					   std::vector<int> &temp_cons_ids)
{
  // we could get clever here and take advantage of the fact that certain
  // quantities are both conserved and primitives (e.g. density and B-field)
  // But, special care needs to be taken. (If we match the primitive density to
  // the conserved density, we will need to modify our definition of primitive
  // density at the half-time step to refer to the temporary density)

  // First get the conserved ids
  

  
  // primitive cell-centered fields:
  //     velocity_x, velocity_y, velocity_z, pressure (thermal)
  // half-timestep cell-centered conserved quantities
  //     temp_density, temp_momentum, temp_energy, temp_B-field
  // half-timestep face-centered longitudinal B-field

  // For each dimension, construct temporary flux fields centered along the
  // faces in the given dimension
  // -- For each dimension, save 1 field for each conserved quantity
  //    (nominally 8 fields)
  // -- We could cheat a little and create a dummy field that is used to
  //    represent the flux along the longitudinal bfield component

  // Allocate the temporary l/r primitive face-centered fields:
  // -- density_l/r, velocity_x_l/r, velocity_y_l/r, velocity_z_l/r,
  //    pressure_l/r bfield_y_l/r, bfield_z_l/r
  // -- We are going to cheat a little and say that the fields are
  //    face-centered along every dimension (so we can reuse the
  //    temporary fields)


  // Allocate the 3 temporary fields for the cell-centered electric fields
  // and the edge-centered Electric fields

  // implement
}

//----------------------------------------------------------------------

void EnzoMethodVlct::deallocate_temp_fields_(Block *block,
					     std::vector<int> &prim_ids,
					     std::vector<int> &priml_ids,
					     std::vector<int> &primr_ids,
					     std::vector<int> &flux_ids,
					     std::vector<int> &efield_ids,
					     int &center_efield_id,
					     std::vector<int> &temp_cons_ids)
{
  // To be filled in
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
