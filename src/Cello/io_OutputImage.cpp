// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputImage.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @brief    Implementation of the OutputImage class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

OutputImage::OutputImage(int index,
			 const Factory * factory,
			 int process_count,
			 int nrows, int ncols) throw ()
  : Output(index,factory),
    data_(),
    op_reduce_(reduce_sum),
    axis_(axis_z),
    nrows_(nrows),
    ncols_(ncols),
    png_(0)

{
  // Override default Output::process_stride_: only root writes
  set_process_stride(process_count);

  map_r_.resize(2);
  map_g_.resize(2);
  map_b_.resize(2);
  map_a_.resize(2);

  map_r_[0] = 0.0;
  map_g_[0] = 0.0;
  map_b_[0] = 0.0;
  map_a_[0] = 1.0;

  map_r_[1] = 1.0;
  map_g_[1] = 1.0;
  map_b_[1] = 1.0;
  map_a_[1] = 1.0;
}

//----------------------------------------------------------------------

OutputImage::~OutputImage() throw ()
{
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void OutputImage::pup (PUP::er &p)
{
  TRACEPUP;

  Output::pup(p);

  p | map_r_;
  p | map_g_;
  p | map_b_;
  p | map_a_;
  p | op_reduce_;
  p | axis_;
  p | nrows_;
  p | ncols_;
  WARNING("OutputImage::pup","skipping data_");
  //  PUParray(p,data_,nrows_*ncols_);
  if (p.isUnpacking()) data_ = 0;
  WARNING("OutputImage::pup","skipping png");
  // p | *png_;
  if (p.isUnpacking()) png_ = 0;
}
#endif

//----------------------------------------------------------------------

void OutputImage::set_colormap
(int n, double * map_r, double * map_g, double * map_b, double * map_a) throw()
{
  // Copy rbg lists

  map_r_.resize(n);
  map_g_.resize(n);
  map_b_.resize(n);
  map_a_.resize(n);

  for (int i=0; i<n; i++) {
    map_r_[i] = map_r[i];
    map_g_[i] = map_g[i];
    map_b_[i] = map_b[i];
    map_a_[i] = map_a ? map_a[i] : 1.0;
  }

}

//======================================================================

void OutputImage::init () throw()
{
  TRACE("OutputImage::init()");
  image_create_();
}

//----------------------------------------------------------------------

void OutputImage::open () throw()
{
  TRACE("OutputImage::open()");
  // Open file if writing a single block
  std::string file_name = expand_file_name_ (&file_name_,&file_args_);

  if (is_writer()) {
    // Create png object
    Monitor::instance()->print ("Output","writing image file %s", 
			      file_name.c_str());
    png_create_(file_name);
  }
}

//----------------------------------------------------------------------

void OutputImage::close () throw()
{
  TRACE("OutputImage::close()");
  if (is_writer()) image_write_();
  image_close_();
  png_close_();
}

//----------------------------------------------------------------------

void OutputImage::finalize () throw()
{
  Output::finalize();
}

//----------------------------------------------------------------------

void OutputImage::write_block
(
 const CommBlock * comm_block,
 const FieldDescr * field_descr
) throw()
// @param comm_block  Block to write
// @param field_descr  Field descriptor
{

  TRACE("OutputImage::write_block()");
  const FieldBlock * field_block = comm_block->block()->field_block();
  
  int nbx,nby,nbz;
  field_block->size(&nbx,&nby,&nbz);

  int ix,iy,iz;

  comm_block->index_forest(&ix,&iy,&iz);

  int ixb0 = ix*nbx;
  int iyb0 = iy*nby;
  int izb0 = iz*nbz;

  // Index of (single) field to write

  it_field_->first();

  int index = it_field_->value();

    // Get ghost depth

  int gx,gy,gz;
  field_descr->ghosts(index,&gx,&gy,&gz);

  // Get array dimensions

  int ndx,ndy,ndz;
  ndx = nbx + 2*gx;
  ndy = nby + 2*gy;
  ndz = nbz + 2*gz;

  // Add block contribution to image

  const char * field_unknowns = field_block->field_unknowns(field_descr,index);

  switch (field_descr->precision(index)) {

  case precision_single:

    image_reduce_ (((float *) field_unknowns),
		   ndx,ndy,ndz, 
		   nbx,nby,nbz,
		   ixb0,iyb0,izb0);
    break;

  case precision_double:

    image_reduce_ (((double *) field_unknowns),
		   ndx,ndy,ndz, 
		   nbx,nby,nbz,
		   ixb0,iyb0,izb0);
    break;

  default:

    WARNING1 ("OutputImage::write_block()",
	     "Unsupported Field precision %d",
	      field_descr->precision(index));
    break;

  }

}

//----------------------------------------------------------------------

void OutputImage::write_field_block
(
 const FieldBlock * Fieldblock,  
 const FieldDescr * field_descr,
 int field_index) throw()
{
  WARNING("OutputImage::write_field_block",
	  "This function should not be called");
}

//----------------------------------------------------------------------

void OutputImage::prepare_remote (int * n, char ** buffer) throw()
{
  TRACE("OutputImage::prepare_remote()");
  DEBUG("prepare_remote");

  int size = 0;
  int nx = nrows_;
  int ny = ncols_;

  // Determine buffer size

  size += nx*ny*sizeof(double);
  size += 2*sizeof(int);
  (*n) = size;

  // Allocate buffer (deallocated in cleanup_remote())
  (*buffer) = new char [ size ];

  union {
    char   * c;
    double * d;
    int    * i;
  } p ;

  p.c = (*buffer);

  *p.i++ = nx;
  *p.i++ = ny;

  for (int k=0; k<nx*ny; k++) *p.d++ = data_[k];
  
}

//----------------------------------------------------------------------

void OutputImage::update_remote  ( int n, char * buffer) throw()
{
  TRACE("OutputImage::update_remote()");
  DEBUG("update_remote");

  union {
    char   * c;
    double * d;
    int    * i;
  } p ;

  p.c = buffer;

  int nx = *p.i++;
  int ny = *p.i++;

  for (int k=0; k<nx*ny; k++) data_[k] += *p.d++;

}

//----------------------------------------------------------------------

void OutputImage::cleanup_remote  (int * n, char ** buffer) throw()
{
  TRACE("OutputImage::cleanup_remote()");
  DEBUG("cleanup_remote");
  delete [] (*buffer);
  (*buffer) = NULL;
}


//======================================================================

void OutputImage::png_create_ (std::string filename) throw()
{
  if (is_writer()) {
    const char * file_name = strdup(filename.c_str());
    png_ = new pngwriter(nrows_, ncols_,0,file_name);
    free ((void *)file_name);
  }
}

//----------------------------------------------------------------------

void OutputImage::png_close_ () throw()
{
  if (is_writer()) {
    png_->close();
    delete png_;
    png_ = 0;
  }
}

//----------------------------------------------------------------------

void OutputImage::image_create_ () throw()
{
  DEBUG("image_create_");

  ASSERT("OutputImage::image_create_",
	 "image_ already created",
	 data_ == NULL);

  data_ = new double [nrows_*ncols_];
  TRACE1("new data_ = %p",data_);

  const double min = std::numeric_limits<double>::max();
  const double max = std::numeric_limits<double>::min();

  double value0;

  switch (op_reduce_) {
  case reduce_min: 
    value0 = max;
    break;
  case reduce_max: 
    value0 = min;
    break;
  case reduce_avg: 
  case reduce_sum: 
  default:         
    value0 = 0; 
    break;
  }

  for (int i=0; i<nrows_*ncols_; i++) data_[i] = value0;

}

//----------------------------------------------------------------------

void OutputImage::image_write_ (double min, double max) throw()
{

  DEBUG("image_write");
  // error check min <= max

  ASSERT2("OutputImage::image_write_",
	 "min %g is greater than max %g",
	 min,max, (min <= max));

  // simplified variable names

  int mx = nrows_;
  int my = ncols_;
  int m  = mx*my;

  // Scale by data if min == max (default)

  if (min == max) {
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::min();
  }

  // Ensure min and max fully enclose data
  for (int i=0; i<m; i++) {
    min = MIN(min,data_[i]);
    max = MAX(max,data_[i]);
  }

  TRACE1("image_write_() data_ = %p",data_);
  size_t n = map_r_.size();

  // loop over pixels (ix,iy)

  for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = data_[i];

      double r=0.0,g=0.0,b=0.0,a=0.0;

      if (min <= value && value <= max) {

	// map v to lower colormap index
	size_t k = size_t((n - 1)*(value - min) / (max-min));

	// prevent k == map_.size()-1, which happens if value == max

	if (k > n - 2) k = n-2;

	// linear interpolate colormap values
	double lo = min +  k   *(max-min)/(n-1);
	double hi = min + (k+1)*(max-min)/(n-1);

	// should be in bounds, but force if not due to rounding error
	if (value < lo) value = lo;
	if (value > hi) value = hi;

	// interpolate colormap

	double ratio = (value - lo) / (hi-lo);

	r = (1-ratio)*map_r_[k] + ratio*map_r_[k+1];
	g = (1-ratio)*map_g_[k] + ratio*map_g_[k+1];
	b = (1-ratio)*map_b_[k] + ratio*map_b_[k+1];
	a = (1-ratio)*map_a_[k] + ratio*map_a_[k+1];
      }

      // Plot pixel
      png_->plot(ix+1, iy+1, 0.0, 0.0, 0.0);
      png_->plot_blend(ix+1, iy+1, a, r,g,b);
    }
  }      

}

//----------------------------------------------------------------------

void OutputImage::image_close_ () throw()
{
  ASSERT("OutputImage::image_create_",
	 "image_ already created",
	 data_ != NULL);
  TRACE1("delete data_ = %p",data_);
  delete [] data_;
  data_ = 0;
}

//----------------------------------------------------------------------

template<class T>
void OutputImage::image_reduce_
(T * array, 
 int ndx, int ndy, int ndz,
 int nx,  int ny,  int nz,
 int nx0, int ny0, int nz0) throw()
{
  // Array multipliers

  int nd3[3] = {1, ndx, ndx*ndy}; 

  // Array size

  int n[3]  = {nx,  ny,  nz};
  int n0[3] = {nx0, ny0, nz0};

  // Remap array axes to image axes axis_x,axis_y

  int axis_x = (axis_+1) % 3;  // image x axis
  int axis_y = (axis_+2) % 3;  // image y-axis

  // Array size permuted to match image

  int npx = n[axis_x];
  int npy = n[axis_y];
  int npz = n[axis_];

  // Array start permuted to match image

  int npx0 = n0[axis_x];
  int npy0 = n0[axis_y];

  // Loop over array subsection

  // image x-axis

  TRACE1("image_reduce_() data_ = %p",data_);
  for (int iax=0; iax<npx; iax++) {

    int iix = npx0 + iax;

    // image y-axis

    for (int iay=0; iay<npy; iay++) {
      
      int iiy = npy0 + iay;

      int index_image = iix + nrows_*iiy;

      if ( ! ( (iix < nrows_) &&
	       (iiy < ncols_)) ) {
	ERROR5 ("OutputImage::image_reduce_",
		"Invalid Access axis %d index(%d %d)  image(%d %d)\n",
		axis_, iix, iiy, nrows_,ncols_);
      }

      double & pixel_value = data_ [index_image];

      // reduce along axis
      for (int iz=0; iz<npz; iz++) {
	
	int index_array = 
	  nd3[axis_x]*iax + 
	  nd3[axis_y]*iay + 
	  nd3[axis_]*iz;

	// reduce along z axis

	switch (op_reduce_) {
	case reduce_min: 
	  pixel_value = MIN(array[index_array],(T)(pixel_value)); 
	  break;
	case reduce_max: 
	  pixel_value = MAX(array[index_array],(T)(pixel_value)); 
	  break;
	case reduce_avg: 
	case reduce_sum: 
	  // @@@@@@@@@@@@@@@@@@@@@@@@@
	  // BUG: segfault here when Charm++ load balancing
	  // @@@@@@@@@@@@@@@@@@@@@@@@@
	  pixel_value += array[index_array]; break;
	default:
	  break;
	}
      }

      if (op_reduce_ == reduce_avg) pixel_value /= npz;

    }

  }
}

