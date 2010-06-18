// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      monitor_Monitor.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Feb 21 12:45:56 PST 2009
/// @todo      simplify image call
/// @brief     Routines for simple output of text, plots, and graphs

#include "cello.hpp"
#include "pngwriter.h"
#include "parallel.hpp"
#include "monitor.hpp" 
#include "error.hpp" 

Monitor * Monitor::instance_ = 0; // (singleton design pattern)

void Monitor::header ()
{
  //    print ("");
  //    print ("     The Laboratory for Computational Astrophysics proudly presents:");
  print ("");
  print ("    =================================================================");
  print ("");
  print ("    oooooooooooo                                          ooooo ooooo ");
  print ("    `888'     `8                                          `888' `888' ");
  print ("     888         ooo. .oo.     oooooooo  .ooooo.           888   888  ");
  print ("     888oooo8    `888P\"Y88b   d'\"\"7d8P  d88' `88b          888   888  ");
  print ("     888    \"     888   888     .d8P'   888   888 8888888  888   888  ");
  print ("     888       o  888   888   .d8P'  .P 888   888          888   888  ");
  print ("    o888ooooood8 o888o o888o d8888888P  `Y8bod8P'         o888o o888o");
  print ("");
  print ("    =================================================================");
  print ("              E N Z O : T H E   N E X T  G E N E R A T I O N");
  print ("    =================================================================");
  print ("");
}

//----------------------------------------------------------------------

void Monitor::image
(std::string name, 
 void * array, 
 precision_type precision,
 int nx, int ny, int nz,
 int nx0, int ny0, int nz0,
 int nx1, int ny1, int nz1,
 int axis, reduce_type op_reduce,
 double min, double max, 
 const double * map, 
 int map_length)
/**
 *********************************************************************
 *
 * @param  name         File name
 * @param  array        Array of values to plot
 * @param  precision    Precision of array elements
 * @param  nx,ny,nz     Size of the array
 * @param  nx0,ny0,nz0  Lower corner of the sub-array
 * @param  nx1,ny1,nz1  Upper bound on the sub-array
 * @param  axis         Which axis to reduce
 * @param  op_reduce    Reduction operator
 * @param  min,max      Bounds for color map values
 * @param  map          Color map [r0, g0, b0, r1, g1, b1, ...]
 * @param  map_length   Length of color map / 3
 *
 * Plot an array as a png file
 *
 *********************************************************************
 */
{

  if (! active_) return;

  // Array size
  int n3[3] = {1, nx, nx*ny}; // Array multipliers

  // Sub-array size
  int m3[3] = {nx1-nx0, ny1-ny0, nz1-nz0};  // sub-array size
  int m0[3] = {nx0, ny0, nz0};              // sub-array start

  // Map array axes to image axes iax,iay
  int iax = (axis+1) % 3;  // image x axis
  int iay = (axis+2) % 3;  // image y-axis
  int iaz = axis;          // reduction axis

  // Array range
  int mx = m3[iax];
  int my = m3[iay];
  int mz = m3[iaz];

  // Array start
  int mx0 = m0[iax];
  int my0 = m0[iay];
  int mz0 = m0[iaz];

  double * image = new double [mx*my];

  // Loop over array subsection

  // image x-axis

  int imax = 0;
  int imin = n3[0]*n3[1]*n3[2];
  for (int jx=0; jx<mx; jx++) {

    int ix = jx + mx0;

    // image y-axis

    for (int jy=0; jy<my; jy++) {
      int iy = jy + my0;
      
      int i = n3[iax]*ix + n3[iay]*iy;
      int j = jx + mx*jy;

      // reduction axis

      long double value = cello::precision_array_value(array,i,precision);

      // initialize reduction
      switch (op_reduce) {
      case reduce_min: image[j] = value; break;
      case reduce_max: image[j] = value; break;
      case reduce_avg: image[j] = 0; break;
      case reduce_sum: image[j] = 0; break;
      default:         image[j] = 0; break;
      }

      for (int jz=0; jz<mz; jz++) {
	int iz = jz + mz0;
	i = n3[iax]*ix + n3[iay]*iy + n3[iaz]*iz;
	long double value = cello::precision_array_value(array,i,precision);
	// reduce along iaz axis
	switch (op_reduce) {
	case reduce_min: image[j] = MIN(value,image[j]); break;
	case reduce_max: image[j] = MAX(value,image[j]); break;
	case reduce_avg: image[j] += value; break;
	case reduce_sum: image[j] += value; break;
	default:         break;
	}
      }
      if (op_reduce == reduce_avg) image[j] /= m3[iaz];
    }
  }

  // Adjust min and max bounds if needed

  double rmin = image[0];
  double rmax = image[0];
  for (int i=0; i<mx*my; i++) {
    if (min > image[i]) min = image[i];
    if (max < image[i]) max = image[i];
    rmin = MIN(rmin,image[i]);
    rmax = MAX(rmax,image[i]);
  }

  pngwriter png (mx,my,0,name.c_str());

  // loop over pixels (jx,jy)
  for (int jx = 0; jx<mx; jx++) {
    for (int jy = 0; jy<my; jy++) {

      int j = jx + mx*jy;
      double v = image[j];

      // map v to lower colormap index
      int index = (map_length-1)*(v - min) / (max-min);

      // prevent index == map_length-1, which happens if v == max
      if (index > map_length - 2) index = map_length-2;

      // linear interpolate colormap values
      double lo = min +  index   *(max-min)/(map_length-1);
      double hi = min + (index+1)*(max-min)/(map_length-1);

      // should be in bounds, but force if not due to rounding error
      if (v < lo) v = lo;
      if (v > hi) v = hi;
      if ( ! (lo <= v && v <= hi)) {
	char buffer [ ERROR_MESSAGE_LENGTH ];
	sprintf (buffer,"v = %g is out of range [%g,%g]",v,lo,hi);
	ERROR_MESSAGE("Montor::image",buffer);
      }

      double ratio = (v - lo) / (hi-lo);
      double r = (1-ratio)*map[3*index+0] + ratio*map[3*(index+1)+0];
      double g = (1-ratio)*map[3*index+1] + ratio*map[3*(index+1)+1];
      double b = (1-ratio)*map[3*index+2] + ratio*map[3*(index+1)+2];
      png.plot(jx+1,jy+1,r,g,b);
    }
  }      

  png.close();

  delete [] image;

}

