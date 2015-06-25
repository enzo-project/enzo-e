// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Matrix.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     yyyy-mm-dd
/// @brief    

#include "compute.hpp"

//----------------------------------------------------------------------

void Matrix::residual (int ir, int ib, int ix, Block * block) throw()
{
  if (! block->is_leaf()) return;

  matvec(ir,ix,block);

  Field field = block->data()->field();

  void * B = field.values(ib);
  void * R = field.values(ir);

  int mx,my,mz;
  field.dimensions(0,&mx,&my,&mz);
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  int precision = field.precision(0);

  if      (precision == precision_single)    
    residual_((float *)(R), (float *)(B),
	      mx,my,mz,nx,ny,nz,gx,gy,gz);
  else if (precision == precision_double)    
    residual_((double *)(R), (double *)(B),
	      mx,my,mz,nx,ny,nz,gx,gy,gz);
  else if (precision == precision_quadruple) 
    residual_((long double *)(R), (long double *)(B),
	      mx,my,mz,nx,ny,nz,gx,gy,gz);
  else 
    ERROR1("Matrix::residual()", "precision %d not recognized", precision);
}

//----------------------------------------------------------------------

template <class T>
void Matrix::residual_ (T * r, T * b,
			int mx, int my, int mz,
			int nx, int ny, int nz,
			int gx, int gy, int gz) throw()
{
  for (int iz=gz; iz<nz+gz; iz++) {
    for (int iy=gy; iy<ny+gy; iy++) {
      for (int ix=gx; ix<nx+gx; ix++) {

	const int i=ix + mx*(iy + mz*iz);

	r[i] = b[i] - r[i];

      }
    }
  }
}
//======================================================================

