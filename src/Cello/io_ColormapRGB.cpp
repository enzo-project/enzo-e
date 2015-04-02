// See LICENSE_CELLO file for license and copyright information

/// @file     io_ColormapRGB.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-06-12 14:42:37
/// @brief    Implementation of an RGB color map

#include "io.hpp"

//----------------------------------------------------------------------

ColormapRGB::ColormapRGB() throw ()
  : Colormap(),
    rgb_()
{
}

//----------------------------------------------------------------------

ColormapRGB::ColormapRGB(const ColormapRGB & colormap) throw ()
  : Colormap(colormap)
/// @param     colormap  Object being copied
{
  rgb_ = colormap.rgb_;
}

//----------------------------------------------------------------------

ColormapRGB & ColormapRGB::operator= (const ColormapRGB & colormap) throw ()
/// @param     colormap  Source object of the assignment
/// @return    The target assigned object
{
  Colormap::operator= (colormap);

  rgb_ = colormap.rgb_;

  return *this;
}

//----------------------------------------------------------------------

void ColormapRGB::pup (PUP::er &p)
{
  TRACEPUP;
  Colormap::pup(p);
  // NOTE: change this function whenever attributes change
  p | rgb_;
}

//----------------------------------------------------------------------

ColormapRGB::~ColormapRGB() throw ()
{
}

//======================================================================

void ColormapRGB::apply (double * r, double * g, double * b,
			 int ndx, int ndy, int ndz,
			 int nx,  int ny,  int nz,
			 float * array)
{
  int n = rgb_.size() / 3;

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int iv = ix+ndx*(iy + ndy*iz);
	int ik = ix+nx*(iy + ny*iz);
	double v = array[iv];
	// map v to lower colormap index
	int k =  (n - 1)*(v - min_) / (max_-min_);
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

void ColormapRGB::apply (double * r, double * g, double * b,
			 int ndx, int ndy, int ndz,
			 int nx,  int ny,  int nz,
			 double * array)
{
  int n = rgb_.size() / 3;

  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      for (int iz=0; iz<nz; iz++) {
	int iv = ix+ndx*(iy + ndy*iz);
	int ik = ix+nx*(iy + ny*iz);
	double v = array[iv];
	// map v to lower colormap index
	int k =  (n - 1)*(v - min_) / (max_-min_);
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
