// $Id: monitor.hpp 1261 2010-03-03 00:14:11Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     monitor_Monitor.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-05
/// @brief    Declaration of the Monitor class

#ifndef MONITOR_MONITOR_HPP
#define MONITOR_MONITOR_HPP

#include "parallel.hpp"
#include "pngwriter.h"
#include "error.hpp"

/// @enum     reduce_type
/// @brief    Reduction operator, used for image projections
enum reduce_type {
  reduce_unknown, /// Unknown reduction
  reduce_min,     /// Minimal value along the axis
  reduce_max,     /// Maximal value along the axis
  reduce_avg,     /// Average value along the axis
  reduce_sum      /// Sum of values along the axis
};

class Monitor {

  /// @class    Monitor
  /// @ingroup  Monitor
  /// @todo     Make calling image() easier
  /// @brief    Functions for user monitoring of the execution status
  ///
  /// The Monitor component is used to communicate information about
  /// the running simulation to the user. Information can be output in
  /// several forms, including text files, HTML files, plots, or other
  /// (generally small) image files. Information is assumed to be from
  /// a correctly-running simulation: anomalous errors or warnings are
  /// output by the Error component. It is assumed that stdout is not
  /// used for monitor output, except for possibly displaying header
  /// text with basic information about Cello and the simulation being
  /// run.

public: // interface

  /// Get single instance of the Monitor object
  static Monitor * instance() throw ()
  { 
    // Delayed creation since Parallel must be initialized
    if (Monitor::instance_ == 0) {
      if (Parallel::instance()->is_initialized()) {
	instance_ = new Monitor(Parallel::instance());
      } else {
	ERROR_MESSAGE("Monitor::instance","Monitor::instance() called before Parallel::initialize()");
      }
    }
    return instance_;
  };

  /// Print the Cello header 
  void header ();

  /// Print a message to stdout
  void print (std::string message, FILE * fp = stdout)
  {
    if (active_) fprintf (fp,"%s %6.1f %s\n",
			  Parallel::instance()->name().c_str(),
			  timer_.value(),message.c_str());
  };

  /// Generate a PNG image of an array
  template<class T>
  void image (std::string name, 
	      T * array,
	      int nx,  int ny,  int nz,   // Array dimensions
	      int nx0, int ny0, int nz0,  // lower inclusive subarray indices
	      int nx1, int ny1, int nz1,  // upper exclusive subarray indices
	      int         axis,           // Axis along which to project
	      reduce_type op_reduce,      // Reduction operation along axis
	      double min, double max,     // Limits for color map
	      const double * color_map,   // color map [r0 g0 b0 r1 g1 b1 ...]
	      int            color_length // length of color map / 3
	      );
  
private: // functions

  /// Initialize the Monitor object (singleton design pattern)
  Monitor(Parallel * parallel) 
    : parallel_(parallel),
      active_(parallel->is_root())
  {  
    timer_.start(); 
  }

private: // attributes

  Parallel * parallel_; // Parallel object, used for is_root()
  bool   active_;  // Whether monitoring is activated.  Used for e.g. ip != 0.
  Timer  timer_;   // Timer from Performance
  
  /// Single instance of the Monitor object (singleton design pattern)
  static Monitor * instance_;

};

//----------------------------------------------------------------------

template<class T>
void Monitor::image
(std::string name, 
 T * array, 
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

#endif /* MONITOR_MONITOR_HPP */

