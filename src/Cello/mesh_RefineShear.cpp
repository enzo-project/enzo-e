// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineShear.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 21 16:02:39 PDT 2014
/// @brief    Implementation of RefineShear class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineShear::RefineShear(const FieldDescr * field_descr,
			 double min_refine,
			 double max_coarsen,
			 std::string output) throw ()
  : Refine (min_refine_, max_coarsen, output)
{
}

//----------------------------------------------------------------------

int RefineShear::apply 
(
 Block * block,
 const FieldDescr * field_descr
 ) throw ()
{

  FieldData * field_data = block->data()->field_data();

  bool all_coarsen = true;
  bool any_refine = false;

  int nx,ny,nz;
  field_data->size(&nx,&ny,&nz);

  int rank = nz > 1 ? 3 : (ny > 1 ? 2 : 1);

  Data * data = block->data();
  double xm[3],xp[3];
  data->lower(&xm[0],&xm[1],&xm[2]);
  data->upper(&xp[0],&xp[1],&xp[2]);

  int id_velocity = field_descr->field_id("velocity_x");

  void * velocity_x = 0;
  void * velocity_y = 0;
  void * velocity_z = 0;

  if (rank >= 1) velocity_x = field_data->values("velocity_x");
  if (rank >= 2) velocity_y = field_data->values("velocity_y");
  if (rank >= 3) velocity_z = field_data->values("velocity_z");
   
  int gx,gy,gz;
  field_descr->ghosts(id_velocity, &gx,&gy,&gz);

  int nxd = nx + 2*gx;
  int nyd = ny + 2*gy;
  int nzd = nz + 2*gz;

  if (rank < 1) gx = 0;
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  void * output = initialize_output_(field_data);

  precision_type precision = field_descr->precision(id_velocity);

  switch (precision) {
  case precision_single:
    evaluate_block_((const float*)velocity_x,
		    (const float*)velocity_y,
		    (const float*)velocity_z,
		    (float*)output,
		    nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		    &any_refine,&all_coarsen, rank);
    break;
  case precision_double:
    evaluate_block_((const double*)velocity_x,
		    (const double*)velocity_y,
		    (const double*)velocity_z,
		    (double*)output,
		    nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		    &any_refine,&all_coarsen, rank);
    break;
  default:
    ERROR2("RefineShear::apply",
	   "Unknown precision %d for velocity_x field",
	   precision,0);
    break;
  }

  int refine_result =  any_refine ?  
    adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

  return refine_result;

}

//----------------------------------------------------------------------

template <class T>
void RefineShear::evaluate_block_(const T * u,
				  const T * v,
				  const T * w,
				  T * output,
				  int ndx, int ndy, int ndz,
				  int nx, int ny, int nz,
				  int gx, int gy, int gz,
				  bool *any_refine,
				  bool * all_coarsen, 
				  int rank)
{
  T shear;
  T uy = 0, vz = 0, wx = 0;
  T uz = 0, vx = 0, wy = 0;
  const int kx = 1;
  const int ky = (rank >= 2) ? ndx : 0;
  const int kz = (rank >= 3) ? ndx*ndy : 0;

  // Compute inner-product of shear vector.  Note works for
  // rank = 1, 2, 3 since 

  for (int iz=gz; iz<nz+gz; iz++) {
    for (int iy=gy; iy<ny+gy; iy++) {
      for (int ix=gx; ix<nx+gx; ix++) {
	int i = ix + ndx*(iy + ndy*iz);
	if (rank >= 2) {
	  uy = u[i+ky] - u[i-ky]; uy *= uy;
	  vx = v[i+kx] - v[i-kx]; vx *= vx;
	} 
	if (rank >= 3) {
	  uz = u[i+kz] - u[i-kz]; uz *= uz;
	  vz = v[i+kz] - v[i-kz]; vz *= vz;
	  wx = w[i+kx] - w[i-kx]; wx *= wx;
	  wy = w[i+ky] - w[i-ky]; wy *= wy;
	}
	shear = uy + uz + vx + vz + wx + wy;
	if (shear > min_refine_)  *any_refine  = true;
	if (shear > max_coarsen_) *all_coarsen = false;
	if (output) {
	  if (shear > max_coarsen_) output[i] =  0;
	  if (shear > min_refine_)  output[i] = +1;
	}
      }
    }
  }
}
//======================================================================

