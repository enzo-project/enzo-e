// See License_CELLO file for license and copyright information

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
			 int nx0, int ny0, int nz0,
			 int nxb, int nyb, int nzb,
			 int max_level,
			 std::string image_type,
			 int image_size_x, int image_size_y,
			 std::string image_reduce_type,
			 std::string image_mesh_color,
			 int         image_block_size,
			 int face_rank,
			 bool image_log,
			 bool ghost,
			 bool specify_bounds,
			 double min, double max) throw ()
  : Output(index,factory),
    image_data_(),
    image_mesh_(),
    axis_(axis_z),
    specify_bounds_(specify_bounds),
    min_(min),max_(max),
    nxi_(image_size_x),
    nyi_(image_size_y),
    png_(0),
    image_type_(image_type),
    face_rank_(face_rank),
    image_log_(image_log),
    ghost_(ghost),
    max_level_(max_level)

{
  if      (image_reduce_type=="min") { op_reduce_ = reduce_min; } 
  else if (image_reduce_type=="max") { op_reduce_ = reduce_max; }
  else if (image_reduce_type=="avg") { op_reduce_ = reduce_avg; }
  else if (image_reduce_type=="sum") { op_reduce_ = reduce_sum; }
  else if (image_reduce_type=="set") { op_reduce_ = reduce_set; }
  else {
    ERROR1 ("OutputImage::OutputImage()",
	    "Unrecognized output_image_reduce_type %s",
	    image_reduce_type.c_str());
  }

  if      (image_mesh_color=="level")   mesh_color_type_ = mesh_color_level;
  else if (image_mesh_color=="process") mesh_color_type_ = mesh_color_process;
  //  else if (image_mesh_color=="neighbor") mesh_color_type_ = mesh_color_neighbor;
  else {
    ERROR1 ("OutputImage::OutputImage()",
	    "Unrecognized output_image_mesh_color %s",
	    image_mesh_color.c_str());
  }

  TRACE1 ("OutputImage reduce %d",op_reduce_);

  int nl = image_block_size * (1 << max_level_); // image size factor

  TRACE1 ("image_block_size factor = %d",nl);

  if (ghost_) {
    if (nx0>1) nx0*=2;
    if (ny0>1) ny0*=2;
    if (nz0>1) nz0*=2;
  }
  if (nxi_ == 0) nxi_ = (type_is_mesh()) ? 2*nl * nxb + 1 : nl * nx0;
  if (nyi_ == 0) nyi_ = (type_is_mesh()) ? 2*nl * nyb + 1 : nl * ny0;
  if (nzi_ == 0) nzi_ = (type_is_mesh()) ? 2*nl * nzb + 1 : nl * nz0;

  TRACE2("OutputImage nl,max_level %d %d",nl,max_level);
  
  // Override default Output::process_stride_: only root writes
  set_process_stride(process_count);

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

void OutputImage::pup (PUP::er &p)
{
  TRACEPUP;

  Output::pup(p);

  p | map_r_;
  p | map_g_;
  p | map_b_;
  WARNING("OutputImage::pup","skipping image_data_");
  if (p.isUnpacking()) image_data_ = 0;
  WARNING("OutputImage::pup","skipping image_mesh_");
  if (p.isUnpacking()) image_mesh_ = 0;
  p | op_reduce_;
  p | mesh_color_type_;
  p | axis_;
  p | specify_bounds_;
  p | min_;
  p | max_;
  p | nxi_;
  p | nyi_;
  p | nzi_;
  WARNING("OutputImage::pup","skipping png");
  // p | *png_;
  if (p.isUnpacking()) png_ = 0;
  p | image_type_;
  p | face_rank_;
  p | image_log_;
  p | ghost_;
  p | max_level_;
}

//----------------------------------------------------------------------

void OutputImage::set_colormap
(int n, double * map_r, double * map_g, double * map_b) throw()
{
  // Copy rbg lists

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
  if (is_writer()) image_write_(min_,max_);
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
 const CommBlock *  comm_block,
 const FieldDescr * field_descr
 ) throw()
// @param comm_block  Block to write
// @param field_descr  Field descriptor
{

  TRACE("OutputImage::write_block()");
  const FieldBlock * field_block = comm_block->block()->field_block();

  // FieldBlock size
  int nbx,nby,nbz;
  field_block->size(&nbx,&nby,&nbz);

  // Block forest array size
  int ix,iy,iz;
  comm_block->index_forest(&ix,&iy,&iz);

  // Index of (single) field to write

  it_field_->first();

  int index_field = it_field_->value();

  // Get ghost depth

  int ngx,ngy,ngz;
  field_descr->ghosts(index_field,&ngx,&ngy,&ngz);

  // FieldBlock array dimensions
  int ndx,ndy,ndz;
  ndx = nbx + 2*ngx;
  ndy = nby + 2*ngy;
  ndz = nbz + 2*ngz;

  // add block contribution to image

  const char * field = (ghost_) ? 
    field_block->field_values(index_field) :
    field_block->field_unknowns(field_descr,index_field);

  // pixel extents of box
  int ixm,iym,izm; 
  int ixp,iyp,izp;
  extents_img_ (comm_block,&ixm,&ixp,&iym,&iyp,&izm,&izp);

  double xm,ym,zm;
  double xp,yp,zp;
  comm_block->lower(&xm,&ym,&zm);
  comm_block->upper(&xp,&yp,&zp);
  double v = 1.0;
  if (nbx > 1) v *= (xp-xm);
  if (nby > 1) v *= (yp-ym);
  if (nbz > 1) v *= (zp-zm);
  TRACE4("output-debug %f %f %f  %f",(xp-xm),(yp-ym),zp-zm,v);

  if (type_is_data() && comm_block->is_leaf()) {

    int mx = ghost_ ? ndx : nbx;
    int my = ghost_ ? ndy : nby;
    int mz = ghost_ ? ndz : nbz;
    for (int ix=0; ix<mx; ix++) {
      int jxm = ixm +  ix   *(ixp-ixm)/mx;
      int jxp = ixm + (ix+1)*(ixp-ixm)/mx-1;
      for (int iy=0; iy<my; iy++) {
	int jym = iym +  iy   *(iyp-iym)/my;
	int jyp = iym + (iy+1)*(iyp-iym)/my-1;
	for (int iz=0; iz<mz; iz++) {
	  // int jzm = izm + iz    *(izp-izm)/mz;
	  // int jzp = izm + (iz+1)*(izp-izm)/mz-1;
	  int i=ix + ndx*(iy + ndy*iz);

	  switch (field_descr->precision(index_field)) {

	  case precision_single:
	    reduce_cube_(image_data_,jxm,jxp,jym,jyp, (((float*)field)[i]));
	    break;

	  case precision_double:
	    reduce_cube_(image_data_,jxm,jxp,jym,jyp,(((double *)field)[i]));
	    break;
	  }
	}
      }
    }
  }

  if (type_is_mesh()) {

    double alpha = 1.0;

    // value for mesh
    double value = 0;
    value = mesh_color_(comm_block->level());
					     
    reduce_box_(image_mesh_,ixm,ixp,iym,iyp,value);

    if (comm_block->is_leaf()) { // ) {
      int xm=(ixm+ixp)/2;
      int ym=(iym+iyp)/2;
      if (face_rank_ <= 1) {
	{
	  int if3[3] = {-1,0,0};
	  int face_level = comm_block->face_level(if3);
	  double face_color = mesh_color_(face_level);
	  reduce_cube_(image_mesh_,ixm+1,ixm+2,ym-1,ym+1, face_color);
	}
	{
	  int if3[3] = {1,0,0};
	  int face_level = comm_block->face_level(if3);
	  double face_color = mesh_color_(face_level);
	  reduce_cube_(image_mesh_,ixp-2,ixp-1,ym-1,ym+1, face_color);
	}
	{
	  int if3[3] = {0,-1,0};
	  int face_level = comm_block->face_level(if3);
	  double face_color = mesh_color_(face_level);
	  reduce_cube_(image_mesh_,xm-1,xm+1,iym+1,iym+2, face_color);
	}
	{
	  int if3[3] = {0,1,0};
	  int face_level = comm_block->face_level(if3);
	  double face_color = mesh_color_(face_level);
	  reduce_cube_(image_mesh_,xm-1,xm+1,iyp-2,iyp-1, face_color);
	}
      }
      if (face_rank_ <= 0) {
	{
	  int if3[3] = {-1,-1,0};
	  int face_level = comm_block->face_level(if3);
	  double face_color = mesh_color_(face_level);
	  reduce_cube_(image_mesh_,ixm+1,ixm+2,iym+1,iym+2, face_color);
	}
	{
	  int if3[3] = {1,-1,0};
	  int face_level = comm_block->face_level(if3);
	  double face_color = mesh_color_(face_level);
	  reduce_cube_(image_mesh_,ixp-2,ixp-1,iym+1,iym+2, face_color);
	}
	{
	  int if3[3] = {-1,1,0};
	  int face_level = comm_block->face_level(if3);
	  double face_color = mesh_color_(face_level);
	  reduce_cube_(image_mesh_,ixm+1,ixm+2,iyp-2,iyp-1, face_color);
	}
	{
	  int if3[3] = {1,1,0};
	  int face_level = comm_block->face_level(if3);
	  double face_color = mesh_color_(face_level);
	  reduce_cube_(image_mesh_,ixp-2,ixp-1,iyp-2,iyp-1, face_color);
	}
      }
    }

  }
}

//----------------------------------------------------------------------

void OutputImage::write_field_block
(
 const FieldBlock * field_block,  
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
  int nx = nxi_;
  int ny = nyi_;

  // Determine buffer size

  size += nx*ny*sizeof(double); // image_data_
  size += nx*ny*sizeof(double); // image_mesh_
  size += 2*sizeof(int);        // nxi_, nyi_
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

  for (int k=0; k<nx*ny; k++) *p.d++ = image_data_[k];

  for (int k=0; k<nx*ny; k++) *p.d++ = image_mesh_[k];
  
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

  if (op_reduce_ == reduce_min) {
    for (int k=0; k<nx*ny; k++) image_data_[k] = std::min(image_data_[k],*p.d++);
    for (int k=0; k<nx*ny; k++) image_mesh_[k] = std::min(image_mesh_[k],*p.d++);
  } else if (op_reduce_ == reduce_max) {
    for (int k=0; k<nx*ny; k++) image_data_[k] = std::max(image_data_[k],*p.d++);
    for (int k=0; k<nx*ny; k++) image_mesh_[k] = std::max(image_mesh_[k],*p.d++);
  } else if (op_reduce_ == reduce_sum) {
    for (int k=0; k<nx*ny; k++) image_data_[k] += *p.d++;
    for (int k=0; k<nx*ny; k++) image_mesh_[k] += *p.d++;
  } else if (op_reduce_ == reduce_avg) {
    for (int k=0; k<nx*ny; k++) image_data_[k]  += *p.d++;
    for (int k=0; k<nx*ny; k++) image_mesh_[k]  += *p.d++;
  } else if (op_reduce_ == reduce_set) {
    for (int k=0; k<nx*ny; k++) image_data_[k]  = *p.d++;
    for (int k=0; k<nx*ny; k++) image_mesh_[k]  = *p.d++;
  }

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

double OutputImage::mesh_color_(int level) const
{
  if (mesh_color_type_ == mesh_color_level) {
    return (1.0+level);
  } else if (mesh_color_type_ == mesh_color_process) {
    return (1.0+CkMyPe())/(CkNumPes());
  } else {
    ERROR1 ("OutputImage::mesh_color_()",
	    "Unknown mesh_color_type_ %d",
	    mesh_color_type_);
    return 0;
  }
}

//----------------------------------------------------------------------

void OutputImage::png_create_ (std::string filename) throw()
{
  if (is_writer()) {
    const char * file_name = strdup(filename.c_str());
    png_ = new pngwriter(nxi_, nyi_,0,file_name);
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
	 image_data_ == NULL || image_mesh_ == NULL);

  image_data_  = new double [nxi_*nyi_];
  image_mesh_  = new double [nxi_*nyi_];
  TRACE2("new image_data_ = %p image_mesh_ = %p",image_data_,image_mesh_);

  const double min = std::numeric_limits<double>::max();
  const double max = std::numeric_limits<double>::min();

  double value0;

  switch (op_reduce_) {
  case reduce_min: 
    value0 = min;
    break;
  case reduce_max: 
    value0 = max;
    break;
  case reduce_avg: 
  case reduce_sum: 
  case reduce_set:
  default:         
    value0 = 0; 
    break;
  }

  for (int i=0; i<nxi_*nyi_; i++) image_data_[i] = value0;
  for (int i=0; i<nxi_*nyi_; i++) image_mesh_[i] = value0;

}

//----------------------------------------------------------------------

void OutputImage::image_write_ (double min_color, double max_color) throw()
{
  //  PARALLEL_PRINTF ("line %d min=%f max=%f\n",__LINE__,min,max);
  // simplified variable names

  int mx = nxi_;
  int my = nyi_;
  int m  = mx*my;

  //  if (! specify_bounds_) {

  double min_val = std::numeric_limits<double>::max();
  double max_val = std::numeric_limits<double>::min();

  //    Compute min and max if needed
  //  if (! specify_bounds_) {
  if (image_log_) {
    for (int i=0; i<m; i++) {
      min_val = MIN(min_val,log(data_(i)));
      max_val = MAX(max_val,log(data_(i)));
    }
  } else {
    for (int i=0; i<m; i++) {
      min_val = MIN(min_val,data_(i));
      max_val = MAX(max_val,data_(i));
    }
  }
  //  }

  double min,max;

  if (specify_bounds_) {
    min = min_color;
    max = max_color;
  } else {
    min = min_val;
    max = max_val;
  }

  TRACE2("image_write_() image_data_ = %p image_mesh_ = %p",image_data_,image_mesh_);
  size_t n = map_r_.size();

  // loop over pixels (ix,iy)

  for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = image_log_ ? log(data_(i)) : data_(i);

      double r=0.0,g=0.0,b=0.0,a=0.0;

	if (value < min) value = min;
	if (value > max) value = max;


      if (min <= value && value <= max) {

	// map v to lower colormap index
	size_t k = size_t((n - 1)*(value - min) / (max-min));

	// prevent k == map_.size()-1, which happens if value == max

	if (k > n - 2) k = n-2;

	// linear interpolate colormap values
	double lo = min +  k   *(max-min)/(n-1);
	double hi = min + (k+1)*(max-min)/(n-1);

	// should be in bounds, but force if not due to rounding error

	// if (value < lo) { 
	//   r=1.0; g=1.0; a=1.0;
	// } else if (value > hi) {
	//   r=0.0; g=0.0; a=0.0;
	// } else {
	// interpolate colormap

	double ratio = (value - lo) / (hi-lo);

	r = (1-ratio)*map_r_[k] + ratio*map_r_[k+1];
	g = (1-ratio)*map_g_[k] + ratio*map_g_[k+1];
	b = (1-ratio)*map_b_[k] + ratio*map_b_[k+1];
	//	if (value < 0.0) { r=1.0; g=0.0; b=0.0; }
	// }
	png_->plot      (ix+1, iy+1, r,g,b);
	//	png_->plot_blend(ix+1, iy+1, a, r,g,b);

      } else {
	
	// red if out of bounds
	png_->plot(ix+1, iy+1, 1.0, 0.0, 0.0);

      }

      // Plot pixel
    }
  }      

}

//----------------------------------------------------------------------

double OutputImage::data_(int index) const
{
  if (type_is_mesh() && type_is_data())
    return image_data_[index] + 0.2*image_mesh_[index];
  else if (type_is_data()) 
    return image_data_[index];
  else  if (type_is_mesh()) 
    return image_mesh_[index];
  else {
    ERROR ("OutputImage::data_()",
	   "image_type is neither mesh nor data");
    return 0.0;
  }
}

//----------------------------------------------------------------------

void OutputImage::image_close_ () throw()
{
  ASSERT("OutputImage::image_create_",
	 "image_ already created",
	 image_data_ != NULL || image_mesh_ != NULL);

  TRACE1("delete image_data_ = %p",image_data_);
  delete [] image_data_;
  image_data_ = 0;

  TRACE1("delete image_mesh_ = %p",image_mesh_)
  delete [] image_mesh_;
  image_mesh_ = 0;
}

//----------------------------------------------------------------------

void OutputImage::reduce_point_ 
(double * data, double value, double alpha) throw()
{
  switch (op_reduce_) {
  case reduce_min:
    *data = std::min(*data,value); 
    break;
  case reduce_max:
    *data = std::max(*data,value); 
    break;
  case reduce_avg:
  case reduce_sum:
    *data += alpha*value;
    break;
  case reduce_set:
    *data = alpha*value;
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_(double * data, int ix0, int ix1, int iy0, int iy1, 
			       double value, double alpha)
{
  if (ix1 < ix0) { int t = ix1; ix1 = ix0; ix0 = t; }
  if (iy1 < iy0) { int t = iy1; iy1 = iy0; iy0 = t; }

  int dx = ix1 - ix0;
  int dy = iy1 - iy0;
  double err = 0.0;
  
  if (dx >= dy) {
    double derr = fabs(1.0*dy/dx);
    int iy = iy0;
    for (int ix = ix0; ix<=ix1; ix++) {
      int i = ix + nxi_*iy;
      reduce_point_(&data[i],value,alpha);
      err += derr;
      if (err >= 0.5) {
	++iy;
	err -= 1.0;
      }
    }
  } else {
    double derr = fabs(1.0*dx/dy);
    int ix = ix0;
    for (int iy = iy0; iy<=iy1; iy++) {
      int i = ix + nxi_*iy;
      reduce_point_(&data[i],value,alpha);
      err += derr;
      if (err >= 0.5) {
	++ix;
	err -= 1.0;
      }
    }
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_x_(double * data, int ixm, int ixp, int iy, double value, double alpha)
{
  ASSERT2("OutputImage::reduce_line_x","! (ixm (%d) <= ixp (%d)",ixm,ixp,ixm<=ixp);
  if (ixp < ixm) { int t = ixp; ixp = ixm; ixm = t; }

  for (int ix=ixm; ix<=ixp; ++ix) {
    int i = ix + nxi_*iy;
    reduce_point_(&data[i],value,alpha);
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_y_(double * data, int ix, int iym, int iyp, double value, double alpha)
{
  ASSERT2("OutputImage::reduce_line_y","! (iym (%d) <= iyp (%d)",iym,iyp,iym<=iyp);
  if (iyp < iym) { int t = iyp; iyp = iym; iym = t; }
  for (int iy=iym; iy<=iyp; ++iy) {
    int i = ix + nxi_*iy;
    reduce_point_(&data[i],value,alpha);
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_box_(double * data,int ixm, int ixp, int iym, int iyp, 
			      double value, double alpha)
{
  TRACE6("reduce_box %d %d %d %d %f %f",ixm,ixp,iym,iyp,value,alpha);
  reduce_line_x_(data,ixm,ixp,iym,value,alpha);
  reduce_line_x_(data,ixm,ixp,iyp,value,alpha);
  reduce_line_y_(data,ixm,iym,iyp,value,alpha);
  reduce_line_y_(data,ixp,iym,iyp,value,alpha);
}

//----------------------------------------------------------------------


void OutputImage::reduce_cube_(double * data, int ixm, int ixp, int iym, int iyp, double value, double alpha)
{
  if (ixm > ixp || iym >iyp) return;

  for (int ix=ixm; ix<=ixp; ++ix) {
    for (int iy=iym; iy<=iyp; ++iy) {
      int i = ix + nxi_*iy;
      reduce_point_(&data[i],value,alpha);
    }
  }
}

//----------------------------------------------------------------------

void OutputImage::extents_img_ (const CommBlock * comm_block,
				int *ixm, int *ixp,
				int *iym, int *iyp,
				int *izm, int *izp) const
{

  int mx,my,mz;

  int ix,iy,iz;
  int nx,ny,nz;
  comm_block->index_global(&ix,&iy,&iz,&nx,&ny,&nz);

  if (image_type_ == "mesh") {
    mx = (nxi_ - 1) / nx;
    my = (nyi_ - 1) / ny;
    mz = (nzi_ - 1) / nz;
  } else {
    mx = (nxi_) / nx;
    my = (nyi_) / ny;
    mz = (nzi_) / nz;
  }
  
  (*ixm) = ix * mx;
  (*iym) = iy * my;
  (*izm) = iz * mz;

  (*ixp) = (ix+1)* mx;
  (*iyp) = (iy+1)* my;
  (*izp) = (iz+1)* mz;
  
}
