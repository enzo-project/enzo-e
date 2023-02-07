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

void pngio::write(const std::string& fname, double* values, int width,
                  int height, const std::vector<float> (&colormap)[3],
                  double vmin, double vmax)
{
  ERROR("pngio::write", "NOT IMPLEMENTED YET");
}

bool* pngio::read_as_mask(const std::string& fname, int& width, int& height)
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
