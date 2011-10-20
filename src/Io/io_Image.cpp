// See LICENSE_CELLO file for license and copyright information

/// @file     image_Image.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-26
/// @brief    Implementation of the Image class

#include "io.hpp"

//----------------------------------------------------------------------

Image::Image() throw ()
  : image_(0),
    nix_(0),
    niy_(0)
{
}

//----------------------------------------------------------------------

void Image::set_map
(int n, double * map_r, double * map_g, double * map_b) throw()
{
  map_r_.resize(n);
  map_g_.resize(n);
  map_b_.resize(n);

  for (int i=0; i<n; i++) {
    map_r_[i] = map_r[i];
    map_g_[i] = map_g[i];
    map_b_[i] = map_b[i];
  }
}

//----------------------------------------------------------------------

void Image::create (int mx, int my) throw()
{
  nix_ = mx;
  niy_ = my;

  ASSERT("Image::image_create_",
	 "image_ already created",
	 image_ == NULL);

  image_ = new double [mx*my];

  for (int i=0; i<mx*my; i++) image_[i] = 0.0;
}

//----------------------------------------------------------------------

void Image::write (pngwriter * png, double min, double max) throw()
{

  // simplified variable names

  int mx = nix_;
  int my = niy_;
  int m  = mx*my;

  // Adjust min and max bounds if needed

   for (int i=0; i<m; i++) {
     min = MIN(min,image_[i]);
     max = MAX(max,image_[i]);
   }

   // loop over pixels (ix,iy)

   for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = image_[i];

      double r = 1.0, g = 0, b = 0;

      if (min <= image_[i] && image_[i] <= max) {

	// map v to lower colormap index
	size_t k = (map_r_.size() - 1)*(value - min) / (max-min);

	// prevent k == map_.size()-1, which happens if value == max
	if (k > map_r_.size() - 2) k = map_r_.size()-2;

	// linear interpolate colormap values
	double lo = min +  k   *(max-min)/(map_r_.size()-1);
	double hi = min + (k+1)*(max-min)/(map_r_.size()-1);

	// should be in bounds, but force if not due to rounding error
	if (value < lo) value = lo;
	if (value > hi) value = hi;

	// interpolate colormap

	double ratio = (value - lo) / (hi-lo);

	r = (1-ratio)*map_r_[k] + ratio*map_r_[k+1];
	g = (1-ratio)*map_g_[k] + ratio*map_g_[k+1];
	b = (1-ratio)*map_b_[k] + ratio*map_b_[k+1];
      }

      // Plot pixel, red if out of range
      png->plot(ix+1, iy+1, r,g,b);
    }
  }      

}

//----------------------------------------------------------------------

void Image::close () throw()
{
  ASSERT("Image::image_create_","image_ already created",
	 image_ != NULL);

  delete [] image_;
  image_ = 0;
}

//======================================================================

