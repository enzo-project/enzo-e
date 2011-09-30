// See LICENSE_CELLO file for license and copyright information

/// @file     io_OutputImage.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @brief    Implementation of the OutputImage class

#include "cello.hpp"

#include "io.hpp"

//----------------------------------------------------------------------

OutputImage::OutputImage(Simulation * simulation) throw ()
  : Output(simulation),
    axis_(axis_z),
    data_(),
    nix_(0),
    niy_(0),
    png_(0)

{
  // Override process_stride_: only root writes

  process_stride_ = simulation->group_process()->size();

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

//----------------------------------------------------------------------

OutputImage::~OutputImage() throw ()
{
}

//----------------------------------------------------------------------

void OutputImage::set_colormap
(int n, double * map_r, double * map_g, double * map_b) throw()
{
  map_r_.resize(n);
  map_g_.resize(n);
  map_b_.resize(n);

  for (int i=0; i<n; i++) {
    map_r_[i] = map_r[i];
    map_g_[i] = map_g[i];
    map_b_[i] = map_b[i];
  }
}

//======================================================================

void OutputImage::init () throw()
{
  // create process image and clear it

  int nxm,nym,nzm;
  simulation_->hierarchy()->patch(0)->size (&nxm, &nym, &nzm);

  image_create_(nxm,nym);

}

//----------------------------------------------------------------------

void OutputImage::open () throw()
{
  // Open file if writing a single block
  std::string file_name = expand_file_name();

  if (is_writer()) {
    int nxm,nym,nzm;
    simulation_->hierarchy()->patch(0)->size (&nxm, &nym, &nzm);
    // Create png object
    Monitor::instance()->print ("[Output] opening image file %s", 
			      file_name.c_str());
    png_create_(file_name,nxm,nym);
  }
}

//----------------------------------------------------------------------

void OutputImage::close () throw()
{
  double min=0.0;
  double max=1.0;
  if (is_writer()) image_write_(min,max);
  image_close_();
  png_close_();
}

//----------------------------------------------------------------------

void OutputImage::write_hierarchy
(
 const FieldDescr * field_descr,
 Hierarchy * hierarchy
  ) throw()
{
  ItPatch it_patch (hierarchy);
  while (Patch * patch = ++it_patch) {

    // NO OFFSET: ASSUMES ROOT PATCH
    write_patch (field_descr, patch,  0,0,0);

  }

}

//----------------------------------------------------------------------

void OutputImage::write_patch
(
 const FieldDescr * field_descr,
 Patch * patch,
 int ixp0,
 int iyp0,
 int izp0
 ) throw()
// @param field_descr  Field descriptor
// @param patch        Patch to output
// @param ixp0  offset of the patch relative to the parent patch along x-axis 
// @param iyp0  offset of the patch relative to the parent patch along y-axis
// @param izp0  offset of the patch relative to the parent patch along z-axis
{

  ItBlock it_block (patch);
  while (Block * block = ++it_block) {

    write_block (field_descr, block, ixp0,iyp0,izp0);
  }

}

//----------------------------------------------------------------------

void OutputImage::write_block
(
 const FieldDescr * field_descr,
 Block * block,
 int ixp0,
 int iyp0,
 int izp0
) throw()
// @param field_descr  Field descriptor
// @param patch        Patch to output
// @param ixp0  offset of the parent patch relative to its parent along x-axis 
// @param iyp0  offset of the parent patch relative to its parent along y-axis 
// @param izp0  offset of the parent patch relative to its parent along z-axis 
{

  FieldBlock * field_block = block->field_block();
  
  int nxb,nyb,nzb;
  field_block->size(&nxb,&nyb,&nzb);

  int ix,iy,iz;
  block->index_patch(&ix,&iy,&iz);

  int ixb0 = ixp0 + ix*nxb;
  int iyb0 = iyp0 + iy*nyb;
  int izb0 = izp0 + iz*nzb;

  // Index of (single) field to write

  it_field_->first();

  int index = it_field_->value();

    // Get ghost depth

  int gx,gy,gz;
  field_descr->ghosts(index,&gx,&gy,&gz);

  // Get array dimensions

  int ndx,ndy,ndz;
  ndx = nxb + 2*gx;
  ndy = nyb + 2*gy;
  ndz = nzb + 2*gz;

  // Add block contribution to image

  char * field_unknowns = field_block->field_unknowns(field_descr,index);

  switch (field_descr->precision(index)) {

  case precision_single:

    image_reduce_ (((float *) field_unknowns),
		   ndx,ndy,ndz, 
		   nxb,nyb,nzb,
		   ixb0,iyb0,izb0,
		   axis_, reduce_sum);
    break;

  case precision_double:

    image_reduce_ (((double *) field_unknowns),
		   ndx,ndy,ndz, 
		   nxb,nyb,nzb,
		   ixb0,iyb0,izb0,
		   axis_, reduce_sum);
    break;

  default:

    WARNING1 ("OutputImage::write_block",
	     "Unsupported Field precision %d",
	      field_descr->precision(index));
    break;

  }

}

//----------------------------------------------------------------------

void OutputImage::prepare_remote (int * n, char ** buffer) throw()
{

  int size = 0;
  int nx = nix_;
  int ny = niy_;

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

void OutputImage::cleanup_remote  ( int * n, char ** buffer) throw()
{
  delete [] (*buffer);
  (*buffer) = NULL;
}


//======================================================================

void OutputImage::png_create_ (std::string filename, int mx, int my) throw()
{
  if (is_writer()) {
    png_ = new pngwriter(mx,my,0,filename.c_str());
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

void OutputImage::image_create_ (int mx, int my) throw()
{
  nix_ = mx;
  niy_ = my;

  ASSERT("OutputImage::image_create_",
	 "image_ already created",
	 data_ == NULL);

  data_ = new double [mx*my];

  for (int i=0; i<mx*my; i++) data_[i] = 0.0;
}

//----------------------------------------------------------------------

void OutputImage::image_write_ (double min, double max) throw()
{

  // simplified variable names

  int mx = nix_;
  int my = niy_;
  int m  = mx*my;

  // Adjust min and max bounds if needed

   for (int i=0; i<m; i++) {
     min = MIN(min,data_[i]);
     max = MAX(max,data_[i]);
   }

   // loop over pixels (ix,iy)

   for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = data_[i];

      double r = 1.0, g = 0, b = 0;

      if (min <= value && value <= max) {

	// map v to lower colormap index
	size_t k = (map_r_.size() - 1)*(value - min) / (max-min);

	// prevent k == map_.size()-1, which happens if value == max
	if (k > map_r_.size() - 2) k = map_r_.size()-2;

	// linear interpolate colormap values
	double lo = min +  k   *(max-min)/(map_r_.size()-1);
	double hi = min + (k+1)*(max-min)/(map_r_.size()-1);

	// should be in bounds, but force if not due to rounding error
	if (value < lo) value = lo;
	if (value > hi) value = hi;

	// interpolate colormap

	double ratio = (value - lo) / (hi-lo);

	r = (1-ratio)*map_r_[k] + ratio*map_r_[k+1];
	g = (1-ratio)*map_g_[k] + ratio*map_g_[k+1];
	b = (1-ratio)*map_b_[k] + ratio*map_b_[k+1];
      }

      // Plot pixel, red if out of range
      if (is_writer()) png_->plot(ix+1, iy+1, r,g,b);
    }
  }      

}

//----------------------------------------------------------------------

void OutputImage::image_close_ () throw()
{
  ASSERT("OutputImage::image_create_",
	 "image_ already created",
	 data_ != NULL);
  delete [] data_;
  data_ = 0;
}

