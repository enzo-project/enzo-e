// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputImage.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Mar 14 17:35:56 PDT 2011
/// @brief    [\ref Io] Declaration for the OutputImage component

#ifndef IO_OUTPUT_IMAGE_HPP
#define IO_OUTPUT_IMAGE_HPP

class ItField;

class OutputImage : public Output {

  /// @class    OutputImage
  /// @ingroup  Io
  /// @brief [\ref Io] class for writing images

public: // functions

  /// Create an uninitialized OutputImage object
  OutputImage(Simulation * simulation) throw();

  /// OutputImage destructor: free allocated image data
  virtual ~OutputImage() throw();

public: // virtual functions

  /// Prepare for accumulating block data
  virtual void init () throw();

  /// Open (or create) a file for IO
  virtual void open () throw();

  /// Close file for IO
  virtual void close () throw();

  /// Write hierarchy-related field data
  virtual void write 
  ( const FieldDescr * field_descr,
    Hierarchy * hierarchy) throw();

  /// Write patch-related field data; may be called by write (Hierarchy)
  virtual void write 
  ( const FieldDescr * field_descr,
    Patch * patch,
    int ix0=0, int iy0=0, int iz0=0) throw();

  /// Write block-related field data; may be called by write (Patch)
  virtual void write 
  ( const FieldDescr * field_descr,
    Block * block,
    int ix0=0, int iy0=0, int iz0=0) throw();

  /// Generate a PNG image of an array
  template<class T>
  void image
  ( std::string name, 
    int mx, int my,             // image size
    T * array,
    int nxd, int nyd, int nzd,   // Array dimensions
    int nx,  int ny,  int nz,   // Array size
    int nx0, int ny0, int nz0,  // Array offset into image
    axis_enum   axis,           // Axis along which to project
    reduce_enum op_reduce,      // Reduction operation along axis
    double min, double max     // Limits for color map
    ) throw();

  // Set the image colormap
  void image_set_map
  (int n, double * map_r, double * map_g, double * map_b) throw();

private: // functions

  /// Create the png file object
  void png_create_ (std::string filename,
		    int image_size_x,  int image_size_y) throw();

  /// Delete the png object
  void png_close_() throw();

  /// Create the image data object
  void image_create_ (int image_size_x,  int image_size_y) throw();

  /// Generate PNG image, using given min and max for colormap
  void image_write_ (double min, double max) throw();

  /// Close the image data
  void image_close_ () throw();

   /// Generate a PNG image of an array
   template<class T>
   void image_reduce_
   ( T * array,
     int nxd, int nyd, int nzd,   // Array dimensions
     int nx,  int ny,  int nz,   // Array dimensions
     int nx0, int ny0, int nz0,  // Array offset into image
     axis_enum   axis,           // Axis along which to project
     reduce_enum op_reduce) throw();


protected: // attributes

  /// Color map
  std::vector<double> map_r_;
  std::vector<double> map_g_;
  std::vector<double> map_b_;

  /// Current image
  double * image_;

  /// Current image size
  int image_size_x_;
  int image_size_y_;

  /// Current pngwriter
  pngwriter * png_;

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

   // Remap array axes to image axes iax,iay

   int iax = (axis+1) % 3;  // image x axis
   int iay = (axis+2) % 3;  // image y-axis
   int iaz = axis;          // reduction axis

   // Array size permuted to match image

   int npx = n[iax];
   int npy = n[iay];
   int npz = n[iaz];

   // Array start permuted to match image

   int npx0 = n0[iax];
   int npy0 = n0[iay];

   // Loop over array subsection

   // image x-axis

   TRACE3("OutputImage::image_reduce_ [%d %d] %d",
	  npx,npy,npz);
   double min=std::numeric_limits<double>::max();
   double max=std::numeric_limits<double>::min();

   for (int index_array_x=0; index_array_x<npx; index_array_x++) {

     int index_image_x = npx0 + index_array_x;

     // image y-axis

     for (int index_array_y=0; index_array_y<npy; index_array_y++) {
      
       int index_image_y = npy0 + index_array_y;

       int index_image = index_image_x + image_size_x_*index_image_y;

       if ( ! ( ( index_image_x < image_size_x_) &&
		 (index_image_y < image_size_y_)) ) {
	 ERROR5 ("OutputImage::image_reduce_",
		 "Invalid Access axis %d index(%d %d)  image(%d %d)\n",
		 axis, index_image_x, index_image_y, image_size_x_,image_size_y_);
       }

       double & pixel_value = image_ [index_image];

       // reduction axis

       // initialize reduction
       switch (op_reduce) {
       case reduce_min: 
	 pixel_value = std::numeric_limits<double>::max();
	 break;
       case reduce_max: 
	 pixel_value = std::numeric_limits<double>::min();
	 break;
       case reduce_avg: 
       case reduce_sum: 
       default:         
	 pixel_value = 0; break;
       }

       // reduce along axis
       for (int iz=0; iz<npz; iz++) {
	
	 int index_array = 
	   nd3[iax]*index_array_x + 
	   nd3[iay]*index_array_y + 
	   nd3[iaz]*iz;

	 // reduce along iaz axis

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
       min=MIN(pixel_value,min);
       max=MAX(pixel_value,max);

     }
   }
   TRACE2 ("OutputImage::image_reduce_ min=%f max=%f",
	   min,max);
 }

//======================================================================

template<class T>
void OutputImage::image
(std::string filename, 
 int mx, int my,
 T * array, 
 int nxd,  int nyd,  int nzd,
 int nx,   int ny,   int nz,
 int nx0, int ny0,   int nz0,
 axis_enum axis, reduce_enum op_reduce,
 double min, double max) throw()

/**
*********************************************************************
*
* @param  filename     File name
* @param  mx,my        Size of the image
* @param  array        Array of values to plot
* @param  nxd,nyd,nzd  Dimension of the array
* @param  nx,ny,nz     Size of the array
* @param  nx0,ny0,nz0  Starting index of the array in the image
* @param  axis         Which axis to reduce
* @param  op_reduce    Reduction operator
* @param  min,max      Bounds for color map values
*
* Plot an array as a png file
*
*********************************************************************
*/
{

  // Return if image is degenerate
  if (mx <= 1 || my <= 1) return;

  // n array size
  // m image size

  // k colormap index

  // Open the image


  png_create_ (filename,mx,my);

  image_create_(mx,my);

  image_reduce_(array,
	       nxd,nyd,nzd,
	       nx,ny,nz,
	       nx0,ny0,nz0,
	       axis,op_reduce);

  // write the image

  image_write_(min,max);

  // delete the image
  image_close_();
}

#endif /* IO_OUTPUT_IMAGE_HPP */
