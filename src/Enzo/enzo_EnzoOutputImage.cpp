// $Id: enzo_EnzoOutputImage.cpp 2093 2011-03-12 01:17:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoOutputImage.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @brief    Implementation of the EnzoOutputImage class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoOutputImage::EnzoOutputImage() throw ()
  : Output(),
    image_(0),
    image_size_x_(0),
    image_size_y_(0),
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

//----------------------------------------------------------------------

EnzoOutputImage::~EnzoOutputImage() throw ()
{
}

//======================================================================

void EnzoOutputImage::write
(
 const FieldDescr * field_descr,
 int index,
 Mesh * mesh,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
  ) throw()
{

  if (mesh->dimension() <= 1) {
    WARNING("EnzoOutputImage::write[Mesh]",
	    "EnzoOutputImage only supports 2D and 3D problems");
  }

  // Open file if writing a single block
  if (root_call) {

    // Get file name
    std::string file_prefix = expand_file_name (cycle,time);
    std::string field_name  = field_descr->field_name(index);
    std::string file_name = file_prefix + "-" + field_name + ".png";

    // Monitor output
    Monitor::instance()->print ("[Output] writing mesh image %s", 
				file_name.c_str());

    // Get mesh size
    int nxm,nym,nzm;
    mesh->patch(0)->size (&nxm, &nym, &nzm);
    // Create image 
    image_open_(file_name,nxm,nym);
  }

  ItPatch it_patch (mesh);
  while (Patch * patch = ++it_patch) {

    // Write patch contribution 
    // NO OFFSET: ASSUMES ROOT PATCH
    write (field_descr, index, patch, mesh, cycle,time,false,
	   0,0,0);
  }

  if (root_call) {
    // Complete geterating image if writing a single block
    double min=0.0;
    double max=1.0;
    if (root_call) image_close_(min,max);
  }

}

//----------------------------------------------------------------------

void EnzoOutputImage::write
(
 const FieldDescr * field_descr,
 int index,
 Patch * patch,
 Mesh  * mesh,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
 ) throw()
{
  if (mesh->dimension() <= 1) {
    WARNING("EnzoOutputImage::write[Patch]",
	    "EnzoOutputImage only supports 2D and 3D problems");
  }

  // Open file if writing a single block
  if (root_call) {

    // Get file name
    std::string file_prefix = expand_file_name (cycle,time);
    std::string field_name = field_descr->field_name(index);
    std::string file_name = file_prefix + "-" + field_name + ".png";

    // Monitor output
    Monitor::instance()->print ("[Output] writing patch image %s", 
				file_name.c_str());

    // Get patch size
    int nxp, nyp;
    patch->size (&nxp, &nyp);

    // Create image 
    image_open_(file_name,nxp,nyp);
  }

#ifdef CONFIG_USE_CHARM
  if (patch->blocks_allocated()) {
    CkPrintf ("%s:%d Output blocks in Patch %p\n",__FILE__,__LINE__,patch);
    //    patch->blocks()
  }
    
#else
  ItBlockLocal it_block (patch);
  while (Block * block = ++it_block) {
    // Get block size
    int nxb,nyb,nzb;
    FieldBlock * field_block = block->field_block();
    field_block->size(&nxb,&nyb,&nzb);

    int ix,iy,iz;
    block->index_patch(&ix,&iy,&iz);

    write (field_descr, index, block, patch, mesh, cycle,time,false,
	   ix0+ix*nxb,
	   iy0+iy*nyb,
	   iz0+iz*nzb);
  }
#endif

  if (root_call) {
    // Complete geterating image if writing a single block
    double min=0.0;
    double max=1.0;
    if (root_call) image_close_(min,max);
  }

}

//----------------------------------------------------------------------

void EnzoOutputImage::write
(
 const FieldDescr * field_descr,
 int index,
 Block * block,
 Patch * patch,
 Mesh  * mesh,
 int cycle,
 double time,
 bool root_call,
 int ix0,
 int iy0,
 int iz0
) throw()
{
  FieldBlock * field_block = block->field_block();

  // Get block size
  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  std::string file_prefix = expand_file_name (cycle,time);

  // Get ghost depth
  int gx,gy,gz;
  field_descr->ghosts(index,&gx,&gy,&gz);

  // Get array dimensions
  int ndx,ndy,ndz;
  ndx=nx+2*gx;
  ndy=ny+2*gy;
  ndz=nz+2*gz;

  if (root_call) {
    // Open file if writing a single block
    std::string field_name = field_descr->field_name(index);
    std::string file_name = file_prefix + "-" + field_name + ".png";
    image_open_(file_name,nx,ny);

    // Monitor output
    Monitor::instance()->print ("[Output] writing block image %s", 
				file_name.c_str());

  }

  // Add block contribution to image

  Scalar * field_unknowns = (Scalar *)field_block->field_unknowns(field_descr,index);
    
  image_reduce_ (field_unknowns, 
		 ndx,ndy,ndz, 
		 nx,ny,nz,
		 ix0,iy0,iz0,
		 axis_z, reduce_sum);

  if (root_call) {
    // Complete geterating image if writing a single block
    double min=0.0;
    double max=1.0;
    if (root_call) image_close_(min,max);
  }
}

//======================================================================

void EnzoOutputImage::image_set_map_
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

//----------------------------------------------------------------------

void EnzoOutputImage::image_open_ (std::string filename, int mx, int my) throw()
{
  png_ = new pngwriter(mx,my,0,filename.c_str());

  image_size_x_ = mx;
  image_size_y_ = my;

  image_ = new double [mx*my];

  for (int i=0; i<mx*my; i++) image_[i] = 0.0;
}

//----------------------------------------------------------------------

void EnzoOutputImage::image_close_ (double min, double max) throw()
{

  // simplified variable names

  int mx = image_size_x_;
  int my = image_size_y_;
  int m  = mx*my;

  // Adjust min and max bounds if needed

  for (int i=0; i<m; i++) {
    min = MIN(min,image_[i]);
    max = MAX(max,image_[i]);
  }

  // loop over pixels (ix,iy)

  for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = image_[i];

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

      double r = (1-ratio)*map_r_[k] + ratio*map_r_[k+1];
      double g = (1-ratio)*map_g_[k] + ratio*map_g_[k+1];
      double b = (1-ratio)*map_b_[k] + ratio*map_b_[k+1];

      png_->plot(ix+1, iy+1, r,g,b);
    }
  }      

  png_->close();

  delete [] image_;
  image_ = 0;
}
