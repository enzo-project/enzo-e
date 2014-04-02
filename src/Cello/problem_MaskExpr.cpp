// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MaskExpr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    Implementation of the MaskExpr class

#include "problem.hpp"

//----------------------------------------------------------------------

MaskExpr::MaskExpr
(Param * param) throw()
  : Mask(), param_(param)
{
}

//----------------------------------------------------------------------

void MaskExpr::copy_(const MaskExpr & mask) throw()
{
  param_ = mask.param_;
}

//----------------------------------------------------------------------

bool MaskExpr::evaluate (double t, double x, double y, double z) const
{
  bool value;
  param_->evaluate_logical(1,&value,&x,&y,&z,t);
  return value;
}

//----------------------------------------------------------------------

void MaskExpr::evaluate (bool * mask, double t,
			 int ndx, int nx, double * xv,
			 int ndy, int ny, double * yv,
			 int ndz, int nz, double * zv) const
{

  ASSERT6("MaskExpr::evaluate",
	  "mask dimension (%d %d %d) needs to be at least (%d %d %d)",
	  ndx,ndy,ndz,nx,ny,nz,
	  (ndx >= nx) && (ndy >= ny) && (ndz >= nz));

  double * x = new double [nx*ny*nz];
  double * y = new double [nx*ny*nz];
  double * z = new double [nx*ny*nz];
  bool * mask_temp = new bool [nx*ny*nz];
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	x[i] = xv[ix];
	y[i] = yv[iy];
	z[i] = zv[iz];
	mask_temp[i] = false;
      }
    }
  }

  int n=nx*ny*nz;

  param_->evaluate_logical(n,mask_temp,x,y,z,t);

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	int id=ix + ndx*(iy + ndy*iz);
	mask[id] = mask_temp[i];
      }
    }
  }

  delete [] mask_temp;
  delete [] z;
  delete [] y;
  delete [] x;

}
