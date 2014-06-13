// See LICENSE_CELLO file for license and copyright information

/// @file     test_Colormap.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-06-12
/// @brief    Test program for the ColormapRGB class

#include "main.hpp"
#include "test.hpp"

#include "io.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  const int nx = 10;
  const int ny = 12;
  const int nz = 3;

  const int ndx = nx + 2;
  const int ndy = ny + 2;
  const int ndz = nz + 2;

  const int nd = ndx*ndy*ndz;

  float  * field = new float [nd];

  float field_min = std::numeric_limits<float>::max();
  float field_max = std::numeric_limits<float>::min();
  for (int ix=1; ix<=nx; ix++) {
    for (int iy=1; iy<=ny; iy++) {
      for (int iz=1; iz<=nz; iz++) {
	int i = ix + ndx*(iy + ndy*iz);
	field[i] = ix+iy+iz;
	field_min = std::min(field[i],field_min);
	field_max = std::max(field[i],field_max);
      }
    }
  }


  //--------------------------------------------------

  // Simple greyscale color map

  double array[] = {0.0,0.0,0.0,  1.0,1.0,1.0};
  std::vector<double> rgb(array,array+(sizeof(array)/sizeof(double)));
  ColormapRGB * colormap = new ColormapRGB(rgb);

  unit_class("ColormapRGB");
  unit_func("ColormapRGB");
  unit_assert (colormap != NULL);

  unit_func("set_limit");
  colormap->set_limit(field_min,field_max);
  unit_func("load");

  colormap->load (ndx,ndy,ndz, nx,ny,nz, &field[1+ndx*(1+ndy*1)]);

  int n = nx*ny*nz;
  double *r, *g, *b;
  r = new double [n];
  g = new double [n];
  b = new double [n];

  unit_func("apply");

  colormap->apply (r,g,b);

  unit_assert(r[0]==0.0);
  unit_assert(g[0]==0.0);
  unit_assert(b[0]==0.0);
  unit_assert(r[n-1]==1.0);
  unit_assert(g[n-1]==1.0);
  unit_assert(b[n-1]==1.0);

  // add more checks of interior values
  unit_assert(false);


  //--------------------------------------------------
  // add more complicated tests
  unit_assert(false);
  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

