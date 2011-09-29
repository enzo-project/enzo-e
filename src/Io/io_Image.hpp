// See LICENSE_CELLO file for license and copyright information

#ifndef IO_IMAGE_HPP
#define IO_IMAGE_HPP

/// @file     io_Image.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-09-26
/// @brief    [\ref Io] Declaration of the Image class
///

class Image {

  /// @class    Image
  /// @ingroup  Io
  /// @brief    [\ref Io] 

public: // interface

  /// Constructor
  Image() throw();

  /// Create the image data object
  void create (int image_size_x,  int image_size_y) throw();

  /// Generate PNG image, using given min and max for colormap
  void write (pngwriter * png, double min, double max) throw();

  /// Close the image data
  void close () throw();

  /// Generate a PNG image of an array
  template<class T>
  void image_reduce_
  ( T * array,
    int nxd, int nyd, int nzd,   // Array dimensions
    int nx,  int ny,  int nz,   // Array dimensions
    int nx0, int ny0, int nz0,  // Array offset into image
    axis_enum   axis,           // Axis along which to project
    reduce_enum op_reduce) throw();

private: // functions

  /// Color map
  std::vector<double> map_r_;
  std::vector<double> map_g_;
  std::vector<double> map_b_;

  /// Current image
  double * image_;

  /// Current image columns
  int nix_;

  /// Current image rows
  int niy_;

private: // attributes


};

//----------------------------------------------------------------------

template<class T>
void OutputImage::image_reduce_
(T * array, 
 int nxd, int nyd, int nzd,
 int nx,  int ny,  int nz,
 int nx0, int ny0, int nz0,
 axis_enum   axis, 
 reduce_enum op_reduce) throw()
{
  // Array multipliers

  int nd3[3] = {1, nxd, nxd*nyd}; 

  // Array size

  int n[3]  = {nx,  ny,  nz};
  int n0[3] = {nx0, ny0, nz0};

  // Remap array axes to image axes axis_x,axis_y

  int axis_x = (axis+1) % 3;  // image x axis
  int axis_y = (axis+2) % 3;  // image y-axis

  // Array size permuted to match image

  int npx = n[axis_x];
  int npy = n[axis_y];
  int npz = n[axis];

  // Array start permuted to match image

  int npx0 = n0[axis_x];
  int npy0 = n0[axis_y];

  // Loop over array subsection

  // image x-axis

  for (int iax=0; iax<npx; iax++) {

    int iix = npx0 + iax;

    // image y-axis

    for (int iay=0; iay<npy; iay++) {
      
      int iiy = npy0 + iay;

      int index_image = iix + nix_*iiy;

      if ( ! ( ( iix < nix_) &&
	       (iiy < niy_)) ) {
	ERROR5 ("OutputImage::image_reduce_",
		"Invalid Access axis %d index(%d %d)  image(%d %d)\n",
		axis, iix, iiy, nix_,niy_);
      }

      double & pixel_value = image_ [index_image];

      double min = std::numeric_limits<double>::max();
      double max = std::numeric_limits<double>::min();

      // reduction axis

      // initialize reduction
      switch (op_reduce) {
      case reduce_min: 
	pixel_value = max;
	break;
      case reduce_max: 
	pixel_value = min;
	break;
      case reduce_avg: 
      case reduce_sum: 
      default:         
	pixel_value = 0; break;
      }

      // reduce along axis
      for (int iz=0; iz<npz; iz++) {
	
	int index_array = 
	  nd3[axis_x]*iax + 
	  nd3[axis_y]*iay + 
	  nd3[axis]*iz;

	// reduce along z axis

	switch (op_reduce) {
	case reduce_min: 
	  pixel_value = MIN(array[index_array],(T)(pixel_value)); 
	  break;
	case reduce_max: 
	  pixel_value = MAX(array[index_array],(T)(pixel_value)); 
	  break;
	case reduce_avg: 
	case reduce_sum: 
	  pixel_value += array[index_array]; break;
	default:
	  break;
	}
      }

      if (op_reduce == reduce_avg) pixel_value /= npz;

    }

  }
}

#endif /* IO_IMAGE_HPP */

