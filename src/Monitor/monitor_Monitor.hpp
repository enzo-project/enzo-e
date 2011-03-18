// $Id$
// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     monitor_Monitor.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-10-05
/// @brief    [\ref Monitor] Declaration of the Monitor class
//----------------------------------------------------------------------
#ifndef MONITOR_MONITOR_HPP
#define MONITOR_MONITOR_HPP

//----------------------------------------------------------------------
/// @def    MONITOR_LENGTH 255
/// @brief  Maximum length of monitor text output

#define MONITOR_LENGTH 255
   
//----------------------------------------------------------------------
/// @enum     reduce_enum
/// @brief    Reduction operator, used for image projections

enum reduce_enum {
  reduce_unknown, /// Unknown reduction
  reduce_min,     /// Minimal value along the axis
  reduce_max,     /// Maximal value along the axis
  reduce_avg,     /// Average value along the axis
  reduce_sum      /// Sum of values along the axis
};

//----------------------------------------------------------------------

class Monitor {

  /// @class    Monitor
  /// @ingroup  Monitor
  /// @todo     Make calling image() easier
  /// @brief    [\ref Monitor] User monitoring of simulation execution status
  ///
  /// The Monitor component is used to communicate information about
  /// the running simulation to the user. Information can be output in
  /// several forms, including text files, HTML files, plots, or other
  /// (generally small) image files. Information is assumed to be from
  /// a correctly-running simulation: anomalous errors or warnings are
  /// output by the Error component.

public: // interface

  static Monitor * instance()
  { if ( instance_ == NULL ) 
      instance_ = new Monitor;
    return instance_;
  };

  /// Set whether the monitor is active for text output.  Useful for
  /// parallel.
  void set_active(bool active) { active_ = active; };

  /// Print the Cello header 
  void header () const;

  /// Print a message to stdout
  void print (std::string buffer, ...) const;

  /// Generate a PNG image of an array
  void image_open (std::string filename, int mx,  int my);

   /// Generate a PNG image of an array
   template<class T>
   void image_write (std::string name, 
   	      T * array,
   	      int nx,  int ny,  int nz,   // Array dimensions
   	      int         axis,           // Axis along which to project
   	      reduce_enum op_reduce,      // Reduction operation along axis
   	      double min, double max     // Limits for color map
		     );

   /// Generate a PNG image of an array
   void image_close ();

  /// Generate a PNG image of an array
  template<class T>
  void image (std::string name, 
	      T * array,
	      int nx,  int ny,  int nz,   // Array dimensions
	      int nx0, int ny0, int nz0,  // Array offset into image
	      int         axis,           // Axis along which to project
	      reduce_enum op_reduce,      // Reduction operation along axis
	      double min, double max     // Limits for color map
	      );

  void set_image_map 
  (int n, double * map_r, double * map_g, double * map_b) throw();

private: // functions

  /// Private constructor of the Monitor object [singleton design pattern]
  Monitor() 
    : active_(true),
      image_(0),
      png_(0)
  { 
    map_r_.resize(2);
    map_g_.resize(2);
    map_b_.resize(2);
    map_r_[0] = 0.0;
    map_g_[0] = 0.0;
    map_b_[0] = 0.0;
    map_r_[1] = 1.0;
    map_g_[1] = 1.0;
    map_b_[1] = 1.0;
  }

  /// Private destructor  of the Monitor object [singleton design pattern]
  ~Monitor()
  {
    delete instance_;
    instance_ = 0;
  }

private: // attributes

  bool   active_;  // Whether monitoring is activated.  Used for e.g. ip != 0.

  /// Color map
  std::vector<double> map_r_;
  std::vector<double> map_g_;
  std::vector<double> map_b_;

  /// Current image
  double * image_;

  /// Current pngwriter
  pngwriter * png_;

  /// Single instance of the Monitor object [singleton design pattern]
  static Monitor * instance_;

};

 //----------------------------------------------------------------------

 template<class T>
 void Monitor::image_write
 (std::string name, 
  T * array, 
  int nx, int ny, int nz,
  int axis, reduce_enum op_reduce,
  double min, double max)
 {
 }

 //----------------------------------------------------------------------

//----------------------------------------------------------------------

template<class T>
void Monitor::image
(std::string filename, 
 T * array, 
 int nx,  int ny,  int nz,
 int nx0, int ny0, int nz0,
 int axis, reduce_enum op_reduce,
 double min, double max)

/**
*********************************************************************
*
* @param  filename     File name
* @param  array        Array of values to plot
* @param  nx,ny,nz     Size of the array
* @param  axis         Which axis to reduce
* @param  op_reduce    Reduction operator
* @param  min,max      Bounds for color map values
*
* Plot an array as a png file
*
*********************************************************************
*/
{
  int nx1, ny1, nz1;

  nx1 = nx0 + nx;
  ny1 = ny0 + ny;
  nz1 = nz0 + nz;

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

  image_open(filename,mx,my);

  // Array start
  int mx0 = m0[iax];
  int my0 = m0[iay];
  int mz0 = m0[iaz];

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
      case reduce_min: image_[j] = array[i]; break;
      case reduce_max: image_[j] = array[i]; break;
      case reduce_avg: image_[j] = 0; break;
      case reduce_sum: image_[j] = 0; break;
      default:         image_[j] = 0; break;
      }

      for (int jz=0; jz<mz; jz++) {
	int iz = jz + mz0;
	i = n3[iax]*ix + n3[iay]*iy + n3[iaz]*iz;
	// reduce along iaz axis
	switch (op_reduce) {
	case reduce_min: image_[j] = MIN(array[i],(T)(image_[j])); break;
	case reduce_max: image_[j] = MAX(array[i],(T)(image_[j])); break;
	case reduce_avg: image_[j] += array[i]; break;
	case reduce_sum: image_[j] += array[i]; break;
	default:         break;
	}
      }
      if (op_reduce == reduce_avg) image_[j] /= m3[iaz];
    }
  }

  // Adjust min and max bounds if needed

  for (int i=0; i<mx*my; i++) {
    min = MIN(min,image_[i]);
    max = MAX(max,image_[i]);
  }


    
  // loop over pixels (jx,jy)
  for (int jx = 0; jx<mx; jx++) {
    for (int jy = 0; jy<my; jy++) {

      int j = jx + mx*jy;
      double v = image_[j];

      // map v to lower colormap index
      int index = (map_r_.size()-1)*(v - min) / (max-min);

      // prevent index == map_.size()-1, which happens if v == max
      if (index > map_r_.size() - 2) index = map_r_.size()-2;

      // linear interpolate colormap values
      double lo = min +  index   *(max-min)/(map_r_.size()-1);
      double hi = min + (index+1)*(max-min)/(map_r_.size()-1);

      // should be in bounds, but force if not due to rounding error
      if (v < lo) v = lo;
      if (v > hi) v = hi;

      double ratio = (v - lo) / (hi-lo);
      double r = (1-ratio)*map_r_[index] + ratio*map_r_[index+1];
      double g = (1-ratio)*map_g_[index] + ratio*map_g_[index+1];
      double b = (1-ratio)*map_b_[index] + ratio*map_b_[index+1];
      png_->plot(jx+1,jy+1,r,g,b);
    }
  }      

  image_close();
}

#endif /* MONITOR_MONITOR_HPP */

