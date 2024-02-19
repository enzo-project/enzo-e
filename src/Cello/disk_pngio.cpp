// See LICENSE_CELLO file for license and copyright information

/// @file     disk_pngio.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-02-06
/// @brief    [\ref Disk] Implementation of functions in the pngio namespace

#include "cello.hpp"
#include "disk.hpp"

// we intentionally include "pngwriter.h" or <png.h> here and NOT in any headers
// so that they don't become transitive dependencies
#include "pngwriter.h"

//----------------------------------------------------------------------

void pngio::write(const std::string& fname, double* data, int width,
                  int height, const std::vector<float> (&colormap)[3],
                  ImgTransform transform,
                  const std::array<double,2>* min_max) noexcept
{

  pngwriter png = pngwriter(width, height, 0, fname.c_str());

  ASSERT("pngio::write", "width must be positive", width > 0);
  ASSERT("pngio::write", "height must be positive", height > 0);

  const int mx = width;
  const int my = height;
  const int m = mx*my;

  const size_t n = colormap[0].size();

  // Determine the min/max values that used to map to colors
  double min = std::numeric_limits<double>::max();
  double max = -std::numeric_limits<double>::max();

  if (min_max != nullptr) {
    // min_max directly specifies the min/max values

    min = MIN(min,(*min_max)[0]);
    max = MAX(max,(*min_max)[1]);
    ASSERT("pngio::write", "max must exceed min", max > min);

  } else {

    if (transform == ImgTransform::log) {
      for (int i=0; i<m; i++) {
        min = MIN(min,log(fabs(data[i])));
        max = MAX(max,log(fabs(data[i])));
      }
    } else if (transform == ImgTransform::abs) {
      for (int i=0; i<m; i++) {
        min = MIN(min,fabs(data[i]));
        max = MAX(max,fabs(data[i]));
      }
    } else if (transform == ImgTransform::none) {
      for (int i=0; i<m; i++) {
        min = MIN(min,data[i]);
        max = MAX(max,data[i]);
      }
    } else {
      ERROR("pngio::write", "transform has unknown value");
    }
  }

  // loop over pixels (ix,iy)

  for (int iy = 0; iy<my; iy++) {
    for (int ix = 0; ix<mx; ix++) {

      int i = ix + mx*iy;

      double value = data[i];

      if (transform == ImgTransform::abs) value = fabs(value);
      if (transform == ImgTransform::log) value = log(fabs(value));

      double r=0.0,g=0.0,b=0.0;

      if (value < min) value = min;
      if (value > max) value = max;

      if (min <= value && value <= max) {

	// map v to lower colormap index
	size_t k =  (n - 1)*(value - min) / (max-min);

	// prevent k == colormap_[0].size()-1, which happens if value == max

	if (k > n - 2) k = n-2;

	// linear interpolate colormap values
	double lo = min +  k   *(max-min)/(n-1);
	double hi = min + (k+1)*(max-min)/(n-1);

	double ratio = (value - lo) / (hi-lo);

	r = (1-ratio)*colormap[0][k] + ratio*colormap[0][k+1];
	g = (1-ratio)*colormap[1][k] + ratio*colormap[1][k+1];
	b = (1-ratio)*colormap[2][k] + ratio*colormap[2][k+1];

	png.plot      (ix+1, iy+1, r,g,b);

      } else {

	// red if out of bounds
	png.plot(ix+1, iy+1, 1.0, 0.0, 0.0);

      }

      // Plot pixel
    }
  }

  // close up the png file (the output data gets written to disk all at once)
  png.close();
}

//----------------------------------------------------------------------

bool* pngio::read_as_mask(const std::string& fname,
                          int& width, int& height) noexcept
{
  pngwriter png;

  // Open the PNG file
  png.readfromfile(fname.c_str());

  // Get the PNG file size

  width = png.getwidth();
  height = png.getheight();

  const int nx = width;
  const int ny = height;
  const int n = nx*ny;

  // Allocate and clear the mask
  bool* mask = new bool [n];
  for (int i=0; i<n; i++) mask[i] = false;

  for (int iy=0; iy<ny; iy++) {
    for (int ix=0; ix<nx; ix++) {

      int i = ix + nx*iy;

      int r = png.read(ix+1,iy+1,1);
      int g = png.read(ix+1,iy+1,2);
      int b = png.read(ix+1,iy+1,3);

      mask[i] = (r+g+b > 0);
    }
  }
  png.close();
  return mask;
}
