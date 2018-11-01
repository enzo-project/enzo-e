
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
}

//----------------------------------------------------------------------

void EnzoMethodVlct::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);
}

//----------------------------------------------------------------------

void EnzoMethodVlct::compute ( Block * block) throw()
{

  if (block->is_leaf()) {

    // To be filled in!
  }

  block->compute_done();
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
