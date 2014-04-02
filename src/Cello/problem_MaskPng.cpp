// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MaskPng.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-04-01
/// @brief    Implementation of the MaskPng class

#include "cello.hpp"

#include "problem.hpp"

//----------------------------------------------------------------------

MaskPng::MaskPng
(std::string file_name, 
 double xm, double xp,
 double ym, double yp) throw()
  : Mask(), mask_(0), nx_(0),ny_(0),
    xm_(xm),xp_(xp),
    ym_(ym),yp_(yp)
{
  pngwriter png;

  // Open the PNG file

  png.readfromfile(file_name.c_str());

  // Get the PNG file size

  nx_ = png.getwidth();
  ny_ = png.getheight();

  // Allocate and clear the mask

  int n = nx_*ny_;
  mask_ = new bool [n];
  for (int i=0; i<n; i++) mask_[i] = false;

  for (int iy=0; iy<ny_; iy++) {
    for (int ix=0; ix<nx_; ix++) {

      int i = ix + nx_*iy;

      int r = png.read(ix+1,iy+1,1);
      int g = png.read(ix+1,iy+1,2);
      int b = png.read(ix+1,iy+1,3);

      mask_[i] = (r+g+b > 0);
    }
  }
  png.close();
}

//----------------------------------------------------------------------

void MaskPng::copy_(const MaskPng & mask) throw()
{
  nx_ = mask.nx_;
  ny_ = mask.ny_;
  const int n = nx_*ny_;
  mask_ = new bool [n];
  for (int i=0; i<n; i++) mask_[i] = mask.mask_[i];
  xm_ = mask.xm_;
  xp_ = mask.xp_;
  ym_ = mask.ym_;
  yp_ = mask.yp_;
}

//----------------------------------------------------------------------

bool MaskPng::evaluate (double t, double x, double y, double z) const
{
  int imx = floor(1.0*( x - xm_) / (xp_-xm_) * nx_);
  int imy = floor(1.0*( y - ym_) / (yp_-ym_) * ny_);
  bool valid = ((0 <= imx && imx < nx_ && 0 <= imy && imy < ny_));
  return valid ? mask_[imx + nx_*imy] : false;
}

//----------------------------------------------------------------------

void MaskPng::evaluate (bool * mask, double t,
			 int ndx, int nx, double * xv,
			 int ndy, int ny, double * yv,
			 int ndz, int nz, double * zv) const
{
  for (int ix=0; ix<nx; ix++) {
    double x = xv[ix];
    int ix_png= floor(1.0*( x - xm_) / (xp_-xm_) * nx_);
    for (int iy=0; iy<ny; iy++) {
      double y = yv[iy];
      int iy_png= floor(1.0*( y - ym_) / (yp_-ym_) * ny_);
      int i_png = ix_png + nx_*(iy_png);
      for (int iz=0; iz<nz; iz++) {
	int i_mask = ix + ndx*(iy + ndy*iz);
	mask[i_mask] = mask_[i_png];
      }
    }
  }
}
