// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_RefineShear.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Jul 21 16:02:39 PDT 2014
/// @brief    Implementation of RefineShear class

#include "mesh.hpp"

//----------------------------------------------------------------------

RefineShear::RefineShear(const FieldDescr * field_descr,
			 double shear_min_refine,
			 double shear_max_coarsen) throw ()
  : shear_min_refine_ (shear_min_refine),
    shear_max_coarsen_(shear_max_coarsen)
{
}

//----------------------------------------------------------------------

int RefineShear::apply 
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

  double h3[3];
  Block * block = comm_block->block();
  double xm[3],xp[3];
  block->lower(&xm[0],&xm[1],&xm[2]);
  block->upper(&xp[0],&xp[1],&xp[2]);

  int id_velocity = field_descr->field_id("velocity_x");

  void * velocity_x = 0;
  void * velocity_y = 0;
  void * velocity_z = 0;

  if (rank >= 1) velocity_x = field_block->values("velocity_x");
  if (rank >= 2) velocity_x = field_block->values("velocity_y");
  if (rank >= 3) velocity_x = field_block->values("velocity_z");
   
  precision_type precision = field_descr->precision(id_velocity);

  switch (precision) {
  case precision_single:
    evaluate_block_((float*)velocity_x,
		    (float*)velocity_y,
		    (float*)velocity_z,
		    nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		    &any_refine,&all_coarsen, rank);
    break;
  case precision_double:
    evaluate_block_((float*)velocity_x,
		    (float*)velocity_y,
		    (float*)velocity_z,
		    nxd,nyd,nzd,nx,ny,nz,gx,gy,gz,
		    &any_refine,&all_coarsen, rank);
    break;
  default:
    ERROR2("RefineShear::apply",
	   "Unknown precision %d for velocity_x field",
	   precision,0);
    break;
  }
}

  int refine_result =  any_refine ?  
    adapt_refine : (all_coarsen ? adapt_coarsen : adapt_same) ;

  return refine_result;

}

//----------------------------------------------------------------------

template <class T>
void RefineShear::evaluate_block_(const T * uin,
				  const T * vin,
				  const T * win,
				  int ndx, int ndy, int ndz,
				  int nx, int ny, int nz,
				  int gx, int gy, int gz,
				  bool *any_refine,
				  bool * all_coarsen, 
				  int rank, 
				  double * h3)
{
  T shear;
  T uy = 0, vz = 0, wx = 0;
  T uz = 0, vx = 0, wy = 0;
  const int kx = 1;
  const int ky = (rank >= 2) ? ndx : 0;
  const int kz = (rank >= 3) ? ndx*ndy : 0;

  // Compute inner-product of shear vector.  Note works for
  // rank = 1, 2, 3 since 
  T temp = 0;
  T * u = (uin?) uin : &temp;
  T * v = (vin?) vin : &temp;
  T * w = (win?) win : &temp;

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i = ix + ndx*(iy + ndy*iz);
	uy = u[i+ky] - u[i-ky]; uy *= uy;
	uz = u[i+kz] - u[i-kz]; uz *= uz;
	vx = v[i+kx] - v[i-kx]; vx *= vx;
	vz = v[i+kz] - v[i-kz]; vz *= vz;
	wx = w[i+kx] - w[i-kx]; wx *= wx;
	wy = w[i+ky] - w[i-ky]; wy *= wy;
	shear = uy + uz + vx + vz + wx + wy;
	if (shear > shear_min_refine_)  *any_refine  = true;
	if (shear > shear_max_coarsen_) *all_coarsen = false;
      }
    }
  }
}
//======================================================================

