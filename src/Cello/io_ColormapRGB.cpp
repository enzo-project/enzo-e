// See LICENSE_CELLO file for license and copyright information

/// @file     io_ColormapRGB.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "io.hpp"

//----------------------------------------------------------------------

ColormapRGB::ColormapRGB() throw ()
  : Colormap(),
    rgb_(),
    ndx_(0),ndy_(0),ndz_(0),
    nx_(0),ny_(0),nz_(0)
{
}

//----------------------------------------------------------------------

ColormapRGB::ColormapRGB(const ColormapRGB & colormap_rgb) throw ()
  : Colormap(colormap_rgb)
/// @param     ColormapRGB  Object being copied
{
}

//----------------------------------------------------------------------

ColormapRGB & ColormapRGB::operator= (const ColormapRGB & colormap_rgb) throw ()
/// @param     ColormapRGB  Source object of the assignment
/// @return    The target assigned object
{
  Colormap::operator= (colormap_rgb);

  return *this;
}

//----------------------------------------------------------------------

void ColormapRGB::pup (PUP::er &p)
{
  TRACEPUP;
  Colormap::pup(p);
  // NOTE: change this function whenever attributes change
  p | rgb_;

  p | ndx_;
  p | ndy_;
  p | ndz_;
  p | nx_;
  p | ny_;
  p | nz_;
}

//----------------------------------------------------------------------

ColormapRGB::~ColormapRGB() throw ()
{
}

//======================================================================

void ColormapRGB::load (int ndx, int ndy, int ndz,
			int nx,  int ny,  int nz,
			float * array)
{
  ndx_ = ndx;
  ndy_ = ndy;
  ndz_ = ndz;
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  precision_ = precision_single;
  vf_ = array;
}

//----------------------------------------------------------------------

void ColormapRGB::load (int ndx, int ndy, int ndz,
			int nx,  int ny,  int nz,
			double * array)
{
  ndx_ = ndx;
  ndy_ = ndy;
  ndz_ = ndz;
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  precision_ = precision_double;
  vd_ = array;
}

//----------------------------------------------------------------------

void ColormapRGB::apply (double * r, double * g, double * b)
{
  int n = rgb_.size() / 3;
  bool is_single = (precision_ == precision_single);

  for (int ix=0; ix<nx_; ix++) {
    for (int iy=0; iy<ny_; iy++) {
      for (int iz=0; iz<nz_; iz++) {
	int iv = ix+ndx_*(iy + ndy_*iz);
	int ik = ix+nx_*(iy + ny_*iz);
	double v = is_single ? vf_[iv] : vd_[iv];
	// map v to lower colormap index
	size_t k =  (n - 1)*(v - min_) / (max_-min_);
	// prevent k == map_.size()-1, which happens if value == max
	if (k > n - 2) k = n-2;
	// linear interpolate colormap values
	double lo = min_ +  k   *(max_-min_)/(n-1);
	double hi = min_ + (k+1)*(max_-min_)/(n-1);

	// should be in bounds, but force if not due to rounding error

	double ratio = (v - lo) / (hi-lo);

	r[ik] = ((1-ratio)*rgb_[3*k+0] + ratio*rgb_[3*(k+1)+0]);
	g[ik] = ((1-ratio)*rgb_[3*k+1] + ratio*rgb_[3*(k+1)+1]);
	b[ik] = ((1-ratio)*rgb_[3*k+2] + ratio*rgb_[3*(k+1)+2]);

      }
    }
  }
}
