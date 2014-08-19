// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRefineShock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 21 16:02:39 PDT 2014
/// @brief    Implementation of EnzoRefineShock class

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoRefineShock::EnzoRefineShock(const FieldDescr * field_descr,
				 double pressure_min_refine,
				 double pressure_max_coarsen,
				 double energy_ratio_min_refine,
				 double energy_ratio_max_coarsen,
				 double gamma,
				 std::string output) throw ()
  : Refine (0.0, 0.0, output),
    pressure_min_refine_ (pressure_min_refine),
    pressure_max_coarsen_(pressure_max_coarsen),
    energy_ratio_min_refine_ (energy_ratio_min_refine),
    energy_ratio_max_coarsen_(energy_ratio_max_coarsen),
    gamma_(gamma)
{
}

//----------------------------------------------------------------------

int EnzoRefineShock::apply 
(
 CommBlock * comm_block,
 const FieldDescr * field_descr
 ) throw ()
{

  FieldBlock * field_block = comm_block->block()->field_block();

  bool all_coarsen = true;
  bool any_refine = false;

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  int rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);

  // compute pressure using the EnzoMethodPressure class
  EnzoMethodPressure method_pressure (EnzoBlock::Gamma);
  method_pressure.compute(comm_block);

  Block * block = comm_block->block();
  double xm[3],xp[3];
  block->lower(&xm[0],&xm[1],&xm[2]);
  block->upper(&xp[0],&xp[1],&xp[2]);

  int id_velocity = field_descr->field_id("velocity_x");

  void * v3[3] = {
    field_block->values("velocity_x"),
    field_block->values("velocity_y"),
    field_block->values("velocity_z") };

  void * te = field_block->values("total_energy");
  void * de = field_block->values("density");
  void * p = field_block->values("pressure");
   
  int gx,gy,gz;
  field_descr->ghosts(id_velocity, &gx,&gy,&gz);

  int nxd = nx + 2*gx;
  int nyd = ny + 2*gy;
  int nzd = nz + 2*gz;

  if (rank < 1) gx = 0;
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  precision_type precision = field_descr->precision(id_velocity);

  void * output = initialize_output_(field_block);

  switch (precision) {
  case precision_single:
    evaluate_block_((const float**) v3,
		    (const float*)  te,
		    (const float*)  de,
		    (const float*)  p,
		    (float*) output,
		    nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		    &any_refine,&all_coarsen, rank);
    break;
  case precision_double:
    evaluate_block_((const double**) v3,
		    (const double*)  te,
		    (const double*)  de,
		    (const double*)  p,
		    (double*) output,
		    nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		    &any_refine,&all_coarsen, rank);
    break;
  default:
    ERROR2("EnzoRefineShock::apply",
	   "Unknown precision %d for velocity_x field",
	   precision,0);
    break;
  }

  int refine_result =  any_refine ?  adapt_refine 
    :                 (all_coarsen ? adapt_coarsen 
		       :             adapt_same) ;

  return refine_result;

}

//----------------------------------------------------------------------

template <class T>
void EnzoRefineShock::evaluate_block_
(const T * v3[],
 const T * te,
 const T * de,
 const T * p,
 T * output,
 int ndx, int ndy, int ndz,
 int nx, int ny, int nz,
 int gx, int gy, int gz,
 bool *any_refine,
 bool * all_coarsen, 
 int rank)
{
  const int d3[3] = {1, ndx, ndx*ndy};

  for (int axis=0; axis<rank; axis++) {

    for (int iz=gz; iz<nz+gz; iz++) {
      for (int iy=gy; iy<ny+gy; iy++) {
	for (int ix=gx; ix<nx+gx; ix++) {

	  int i = ix + ndx*(iy + ndy*iz);
	  int id = d3[axis];

	  T dp = fabs    (p[i+id] - p[i-id]) 
	    / (std::min(p[i+id] , p[i-id])) ;

	  T dv = v3[axis][i+id] - v3[axis][i-id];

	  T e = p[i]/(gamma_ - 1.0);

	  T ep = te[i+id]*de[i+id];
	  T e0 = te[i]   *de[i];
	  T em = te[i-id]*de[i-id];

	  T er = e / std::max (std::max(em,e0),ep);

	  bool l_refine = (dv < 0.0) && 
	    (dp > pressure_min_refine_) &&
	    (er > energy_ratio_min_refine_);

	  bool l_same = (dv < 0.0) &&
	    (dp > pressure_max_coarsen_) &&
	    (er > energy_ratio_max_coarsen_);

	  if (l_refine)  *any_refine = true;
	  if (l_same)    *all_coarsen = false;

	  if (output) {
	    if (l_same)   output[i] =  0;
	    if (l_refine) output[i] = +1;
	  }
	}
      }
    }
  }
}
//======================================================================

