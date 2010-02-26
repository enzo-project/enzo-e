/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */
//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      monitor.cpp
 * @brief     Routines for simple output of text, plots, and graphs
 * @author    James Bordner
 * @date      Thu Feb 21 12:45:56 PST 2009
 * @note      
 *
 * DESCRIPTION 
 * 
 *    Routines for simple output of text, plots, and graphs
 *
 * PACKAGES
 *
 *    none
 * 
 * INCLUDES
 *  
 *    pngwriter.h
 *
 * PUBLIC FUNCTIONS
 *  
 *    
 *
 * PRIVATE FUCTIONS
 *  
 *    
 *
 * $Id$
 *
 *********************************************************************
 */

#include <assert.h>

#include "cello.h"
#include "pngwriter.h"

#include "monitor.hpp" 

void Monitor::image
(std::string name, 
 Scalar * array, 
 int nx, int ny, int nz,
 int nx0, int ny0, int nz0,
 int nx1, int ny1, int nz1,
 int axis, enum_reduce op_reduce,
 Scalar min, Scalar max, 
 const double * map, 
 int map_length)
/**
 *********************************************************************
 *
 * @param  name         File name
 * @param  array        Array of values to plot
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
  for (int jx=0; jx<mx; jx++) {

    int ix = jx + mx0;

    // image y-axis

    for (int jy=0; jy<my; jy++) {
      int iy = jy + my0;
      
      int i = n3[iax]*ix + n3[iay]*iy;
      int j = jx + mx*jy;

      // reduction axis

      // initialize reduction
      switch (op_reduce) {
      case reduce_min: image[j] = array[i]; break;
      case reduce_max: image[j] = array[i]; break;
      case reduce_avg: image[j] = 0; break;
      case reduce_sum: image[j] = 0; break;
      default:         image[j] = 0; break;
      }

      for (int jz=0; jz<mz; jz++) {
	int iz = jz + mz0;
	i = n3[iax]*ix + n3[iay]*iy + n3[iaz]*iz;
	// reduce along iaz axis
	switch (op_reduce) {
	case reduce_min: image[j] = MIN(array[i],image[j]); break;
	case reduce_max: image[j] = MAX(array[i],image[j]); break;
	case reduce_avg: image[j] += array[i]; break;
	case reduce_sum: image[j] += array[i]; break;
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

  double h = (max-min) / (map_length-1);

  pngwriter png (mx,my,0,name.c_str());

  for (int jx = 0; jx<mx; jx++) {
    for (int jy = 0; jy<my; jy++) {

      int j = jx + mx*jy;
      double v = image[j];

      // map v to lower colormap index
      int index = (map_length-1)*(v - min) / (max-min);

      // prevent index == map_length-1, which happens if v == max
      if (index > map_length - 2) index = map_length-2;

      // linear interpolate colormap values
      double lo = min + index*(max-min)/(map_length-1);
      double hi = min + (index+1)*(max-min)/(map_length-1);

      // should be in bounds, but force if not due to rounding error
      if (v < lo) v = lo;
      if (v > hi) v = hi;
      assert (lo <= v && v <= hi);

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

