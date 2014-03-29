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
 int index_parameter,
 double time,
 int nx, double xmin, double xmax,
 int ny, double ymin, double ymax,
 int nz, double zmin, double zmax) throw()
  : parameters_(parameters),
    parameter_name_(parameter_name),
    index_parameter_(index_parameter),
    time_(time),
    nx_(nx),xmin_(xmin),xmax_(xmax),
    ny_(ny),ymin_(ymin),ymax_(ymax),
    nz_(nz),zmin_(zmin),zmax_(zmax)
{
}

//----------------------------------------------------------------------

MaskExpr::~MaskExpr() throw ()
{
}

//----------------------------------------------------------------------

MaskExpr::MaskExpr(const MaskExpr & mask) throw ()
/// @param     MaskExpr  Object being copied
{
  copy_(mask);
}

//----------------------------------------------------------------------

MaskExpr & MaskExpr::operator= (const MaskExpr & mask) throw ()
/// @param     MaskExpr  Source object of the assignment
/// @return    The target assigned object
{
  copy_(mask);
  return *this;
}

//----------------------------------------------------------------------

void MaskExpr::copy_(const MaskExpr & mask) throw()
{
  const int n = mask.nx_*mask.ny_*mask.nz_;

  parameters_ = mask.parameters_;
  parameter_name_ = mask.parameter_name_;
  index_parameter_ = mask.index_parameter_;
  time_ = mask.time_;
  nx_ = mask.nx_;  xmin_ = mask.xmin_;  xmax_ = mask.xmax_;
  ny_ = mask.ny_;  ymin_ = mask.ymin_;  zmin_ = mask.zmin_;
  nz_ = mask.nz_;  ymax_ = mask.ymax_;  zmax_ = mask.ymax_;
}

//----------------------------------------------------------------------

bool MaskExpr::evaluate (int ix, int iy, int iz) const
{
  bool value;
  bool deflt = false;
  double x = xmin_ + (ix+0.5)*(xmax_-xmin_)/nx_;
  double y = ymin_ + (iy+0.5)*(ymax_-ymin_)/ny_;
  double z = zmin_ + (iz+0.5)*(zmax_-zmin_)/nz_;
  double t = time_;
  parameters_->list_evaluate_logical
    (index_parameter_,
     parameter_name_,
     1,
     &value,
     &deflt,
     &x,&y,&z,t);

  return value;
}

//----------------------------------------------------------------------

void MaskExpr::evaluate (bool * mask, int ndx, int ndy, int ndz) const
{

  ASSERT6("MaskExpr::evaluate",
	  "mask dimension (%d %d %d) needs to be at least (%d %d %d)",
	  ndx,ndy,ndz,nx_,ny_,nz_,
	  (ndx >= nx_) && (ndy >= ny_) && (ndz >= nz_));

  double * x = new double [nx_*ny_*nz_];
  double * y = new double [nx_*ny_*nz_];
  double * z = new double [nx_*ny_*nz_];
  double t = time_;
  bool * mask_temp = new bool [nx_*ny_*nz_];
  bool * mask_deflt = new bool [nx_*ny_*nz_];
  for (int ix=0; ix<nx_; ix++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int iz=0; iz<nz_; iz++) {
	int i=ix + nx_*(iy + ny_*iz);
	x[i] = xmin_ + (ix+0.5)*(xmax_-xmin_)/nx_;
	y[i] = ymin_ + (iy+0.5)*(ymax_-ymin_)/ny_;
	z[i] = zmin_ + (iz+0.5)*(zmax_-zmin_)/nz_;
	mask_temp[i] = false;
	mask_deflt[i] = false;
      }
    }
  }

  int n=nx_*ny_*nz_;

  parameters_->list_evaluate_logical
    (index_parameter_,
     parameter_name_,
     n,
     mask_temp,
     mask_deflt,
     x,y,z,t);

  for (int ix=0; ix<nx_; ix++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int iz=0; iz<nz_; iz++) {
	int i=ix + nx_*(iy + ny_*iz);
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
