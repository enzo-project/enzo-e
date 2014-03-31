// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MaskExpr.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-03-28
/// @brief    Implementation of the MaskExpr class

#include <cstring>

#include "problem.hpp"

//----------------------------------------------------------------------

MaskExpr::MaskExpr
(Parameters * parameters,
 const std::string parameter_name,
 int index_parameter) throw()
  : parameters_(parameters),
    parameter_name_(parameter_name),
    index_parameter_(index_parameter)
{
}

//----------------------------------------------------------------------

void MaskExpr::copy_(const MaskExpr & mask) throw()
{
  parameters_ = mask.parameters_;
  parameter_name_ = mask.parameter_name_;
  index_parameter_ = mask.index_parameter_;
}

//----------------------------------------------------------------------

bool MaskExpr::evaluate (double t, double x, double y, double z) const
{
  bool value;
  bool deflt = false;
  parameters_->list_evaluate_logical
    (index_parameter_, parameter_name_, 1, &value, &deflt,  &x,&y,&z,t);

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
  bool * mask_deflt = new bool [nx*ny*nz];
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	x[i] = xv[ix];
	y[i] = yv[iy];
	z[i] = zv[iz];
	mask_temp[i] = false;
	mask_deflt[i] = false;
      }
    }
  }

  int n=nx*ny*nz;

  parameters_->list_evaluate_logical
    (index_parameter_, parameter_name_,  n, mask_temp, mask_deflt, x,y,z,t);

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int i=ix + nx*(iy + ny*iz);
	int id=ix + ndx*(iy + ndy*iz);
	mask[id] = mask_temp[i];
      }
    }
  }

  delete [] mask_deflt;
  delete [] mask_temp;
  delete [] z;
  delete [] y;
  delete [] x;

}
