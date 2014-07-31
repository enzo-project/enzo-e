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

  const int nx = 101;
  const int ny = 1;
  const int nz = 1;

  const int n = nx*ny*nz;

  const int gx = 1;
  const int gy = 0;
  const int gz = 0;

  const int ndx = nx + 2*gx;
  const int ndy = ny + 2*gy;
  const int ndz = nz + 2*gz;

  const int nd = ndx*ndy*ndz;

  float  * array = new float [nd];

  float array_min = std::numeric_limits<float>::max();
  float array_max = std::numeric_limits<float>::min();
  for (int ix=gx; ix<nx+gx; ix++) {
    for (int iy=gy; iy<ny+gy; iy++) {
      for (int iz=gz; iz<nz+gz; iz++) {
	int i = ix + ndx*(iy + ndy*iz);
	array[i] = 2.*ix+7.*iy+19.*iz;
	array_min = std::min(array[i],array_min);
	array_max = std::max(array[i],array_max);
      }
    }
  }

  double * r = new double[n];
  double * g = new double[n];
  double * b = new double[n];

  unit_class("ColormapRGB");

  //--------------------------------------------------
  // Simple greyscale color map

  double rgb_2_array[] = {0.0,0.0,0.0,  1.0,1.0,1.0};
  std::vector<double> rgb_2(rgb_2_array,
			  rgb_2_array+(sizeof(rgb_2_array)/sizeof(double)));
  ColormapRGB * colormap_2 = new ColormapRGB(rgb_2);

  unit_func("ColormapRGB");
  unit_assert (colormap_2 != NULL);

  unit_func("set_limit");
  colormap_2->set_limit(array_min,array_max);

  unit_func("apply");

  colormap_2->apply (r,g,b, ndx,ndy,ndz, nx,ny,nz, &array[gx+ndx*(gy+ndy*gz)]);

  unit_assert(r[ 0 ]==0.0);
  unit_assert(g[ 0 ]==0.0);
  unit_assert(b[ 0 ]==0.0);
  unit_assert(r[n-1]==1.0);
  unit_assert(g[n-1]==1.0);
  unit_assert(b[n-1]==1.0);

  unit_func("set_limit");
  
  colormap_2->set_limit(array_min,2*array_max-array_min);

  colormap_2->apply (r,g,b, ndx,ndy,ndz, nx,ny,nz, &array[gx+ndx*(gy+ndy*gz)]);
  unit_assert(r[ 0 ]==0.0);
  unit_assert(g[ 0 ]==0.0);
  unit_assert(b[ 0 ]==0.0);
  unit_assert(r[n-1]==0.5);
  unit_assert(g[n-1]==0.5);
  unit_assert(b[n-1]==0.5);

  //--------------------------------------------------
  // Mixed colormap

  double rgb_3_array[] = {0.0,0.0,0.0,  1.0,0.5,0.25, 0.0, 0.0, 0.0};
  std::vector<double> rgb_3(rgb_3_array,
			    rgb_3_array+(sizeof(rgb_3_array)/sizeof(double)));
  ColormapRGB * colormap_3 = new ColormapRGB(rgb_3);

  unit_func("ColormapRGB");
  unit_assert (colormap_3 != NULL);

  unit_func("set_limit");
  colormap_3->set_limit(array_min,array_max);

  unit_func("apply");

  colormap_3->apply (r,g,b, ndx,ndy,ndz, nx,ny,nz, &array[gx+ndx*(gy+ndy*gz)]);

  const int mx=(nx+1)/2-1;
  const int my=(ny+1)/2-1;
  const int mz=(nz+1)/2-1;
  const int m = mx+ndx*(my + ndy*mz);

  unit_assert(r[ 0 ]==0.0);
  unit_assert(g[ 0 ]==0.0);
  unit_assert(b[ 0 ]==0.0);
  unit_assert(r[ m ]==1.0);
  unit_assert(g[ m ]==0.5);
  unit_assert(b[ m ]==0.25);
  unit_assert(r[n-1]==0.0);
  unit_assert(g[n-1]==0.0);
  unit_assert(b[n-1]==0.0);

  unit_func("set_limit");
  
  colormap_3->set_limit(array_min,2*array_max-array_min);

  colormap_3->apply (r,g,b, ndx,ndy,ndz, nx,ny,nz, &array[gx+ndx*(gy+ndy*gz)]);
  unit_assert(r[ 0 ]==0.0);
  unit_assert(g[ 0 ]==0.0);
  unit_assert(b[ 0 ]==0.0);
  unit_assert(r[ m ]==0.5);
  unit_assert(g[ m ]==0.25);
  unit_assert(b[ m ]==0.125);
  unit_assert(r[n-1]==1.0);
  unit_assert(g[n-1]==0.5);
  unit_assert(b[n-1]==0.25);

  //--------------------------------------------------

  delete [] r;
  delete [] g;
  delete [] b;


  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

