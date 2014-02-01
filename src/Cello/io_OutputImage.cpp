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
    data_(),
    axis_(axis_z),
    specify_bounds_(specify_bounds),
    min_(min),max_(max),
    nxi_(image_size_x),
    nyi_(image_size_y),
    png_(0),
    image_type_(image_type),
    face_rank_(face_rank),
    image_log_(image_log),
    ghost_(ghost)

{
  //  PARALLEL_PRINTF ("line %d min=%f max=%f\n",__LINE__,min,max);
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

  PARALLEL_PRINTF ("image_mesh_color = %s\n",image_mesh_color.c_str());
  if      (image_mesh_color=="level")   mesh_color_ = mesh_color_level;
  else if (image_mesh_color=="process") mesh_color_ = mesh_color_process;
  else if (image_mesh_color=="neighbor") mesh_color_ = mesh_color_neighbor;
  else {
    ERROR1 ("OutputImage::OutputImage()",
	    "Unrecognized output_image_mesh_color %s",
	    image_mesh_color.c_str());
  }

  TRACE1 ("OutputImage reduce %d",op_reduce_);

  int nl = image_block_size * (1 << max_level); // image size factor

  TRACE1 ("image_block_size factor = %d",nl);

  if (ghost_) {
    if (nx0>1) nx0*=2;
    if (ny0>1) ny0*=2;
    if (nz0>1) nz0*=2;
  }
  if (nxi_ == 0) nxi_ = (image_type_ == "mesh") ? 2*nl * nxb + 1 : nl * nx0;
  if (nyi_ == 0) nyi_ = (image_type_ == "mesh") ? 2*nl * nyb + 1 : nl * ny0;
  if (nzi_ == 0) nzi_ = (image_type_ == "mesh") ? 2*nl * nzb + 1 : nl * nz0;

  TRACE2("OutputImage nl,max_level %d %d",nl,max_level);
  
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

void OutputImage::pup (PUP::er &p)
{
  TRACEPUP;

  Output::pup(p);

  p | map_r_;
  p | map_g_;
  p | map_b_;
  p | map_a_;
  WARNING("OutputImage::pup","skipping data_");
  if (p.isUnpacking()) data_ = 0;
  p | op_reduce_;
  p | mesh_color_;
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
}

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
  //  PARALLEL_PRINTF ("line %d min=%f max=%f\n",__LINE__,min_,max_);
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
 const CommBlock * comm_block,
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

  int gx,gy,gz;
  field_descr->ghosts(index_field,&gx,&gy,&gz);

  // FieldBlock array dimensions
  int ndx,ndy,ndz;
  ndx = nbx + 2*gx;
  ndy = nby + 2*gy;
  ndz = nbz + 2*gz;

  // add block contribution to image

  const char * field = (ghost_) ? 
    field_block->field_values(index_field) :
    field_block->field_unknowns(field_descr,index_field);

  double color = 0;
  if (mesh_color_ == mesh_color_level)   color = comm_block->level()+1;
  if (mesh_color_ == mesh_color_process) color = CkMyPe()+1;
  if (mesh_color_ == mesh_color_neighbor) color = comm_block->count_neighbors();

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

  if (image_type_ == "mesh") {

    reduce_box_(ixm,ixp,iym,iyp,color);

    if (comm_block->is_leaf()) {
      int xm=(ixm+ixp)/2;
      int ym=(iym+iyp)/2;
      if (face_rank_ <= 1) {
	{
	  int if3[3] = {-1,0,0};
	  int face_level = comm_block->face_level(if3)+1;
	  reduce_cube_(ixm+1,ixm+2,ym-1,ym+1, face_level);
	}
	{
	  int if3[3] = {1,0,0};
	  int face_level = comm_block->face_level(if3)+1;
	  reduce_cube_(ixp-2,ixp-1,ym-1,ym+1, face_level);
	}
	{
	  int if3[3] = {0,-1,0};
	  int face_level = comm_block->face_level(if3)+1;
	  reduce_cube_(xm-1,xm+1,iym+1,iym+2, face_level);
	}
	{
	  int if3[3] = {0,1,0};
	  int face_level = comm_block->face_level(if3)+1;
	  reduce_cube_(xm-1,xm+1,iyp-2,iyp-1, face_level);
	}
      }
      if (face_rank_ <= 0) {
	{
	  int if3[3] = {-1,-1,0};
	  int face_level = comm_block->face_level(if3)+1;
	  reduce_cube_(ixm+1,ixm+2,iym+1,iym+2, face_level);
	}
	{
	  int if3[3] = {1,-1,0};
	  int face_level = comm_block->face_level(if3)+1;
	  reduce_cube_(ixp-2,ixp-1,iym+1,iym+2, face_level);
	}
	{
	  int if3[3] = {-1,1,0};
	  int face_level = comm_block->face_level(if3)+1;
	  reduce_cube_(ixm+1,ixm+2,iyp-2,iyp-1, face_level);
	}
	{
	  int if3[3] = {1,1,0};
	  int face_level = comm_block->face_level(if3)+1;
	  reduce_cube_(ixp-2,ixp-1,iyp-2,iyp-1, face_level);
	}
      }
    }

  }
    //  } else if (image_type_ == "data") {

    if (! comm_block->is_leaf()) return;
    // for each cell
    TRACE3("OutputImage ixm,ixp,nbx %d %d %d",ixm,ixp,nbx);
    TRACE3("OutputImage iym,iyp,nby %d %d %d",iym,iyp,nby);
    TRACE3("OutputImage izm,izp,nbz %d %d %d",izm,izp,nbz);

    //@@@@@@@@@@
    // DEBUG
    //@@@@@@@@@@
    // float savef=0.0;
    // double saved=0.0;
    // if (field_descr->precision(index_field) == precision_single) {
    //   savef = ((float*)field)[0];
    //   ((float*)field)[0] = 0.0;
    // } else {
    //   saved = ((double*)field)[0];
    //   ((double*)field)[0] = 0.0;
    // }

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
	    reduce_cube_(jxm,jxp,jym,jyp, (((float*)field)[i]));
	    break;

	  case precision_double:
	    reduce_cube_(jxm,jxp,jym,jyp,(((double *)field)[i]));
	    break;
	  }


	}
      }
    }
    // if (field_descr->precision(index_field) == precision_single) {
    //   ((float*)field)[0] = savef;
    // } else {
    //   ((double*)field)[0] = saved;
    // }
    //  }
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
  int nx = nxi_;
  int ny = nyi_;

  // Determine buffer size

  size += nx*ny*sizeof(double); // data_
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

  if (op_reduce_ == reduce_min) {
    for (int k=0; k<nx*ny; k++) data_[k] = std::min(data_[k],*p.d++);
  } else if (op_reduce_ == reduce_max) {
    for (int k=0; k<nx*ny; k++) data_[k] = std::max(data_[k],*p.d++);
  } else if (op_reduce_ == reduce_sum) {
    for (int k=0; k<nx*ny; k++) data_[k] += *p.d++;
  } else if (op_reduce_ == reduce_avg) {
    for (int k=0; k<nx*ny; k++) data_[k]  += *p.d++;
  } else if (op_reduce_ == reduce_set) {
    for (int k=0; k<nx*ny; k++) data_[k]  = *p.d++;

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
	 data_ == NULL);

  data_  = new double [nxi_*nyi_];
  TRACE1("new data_ = %p",data_);

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

  for (int i=0; i<nxi_*nyi_; i++) data_[i] = value0;

}

//----------------------------------------------------------------------

void OutputImage::image_write_ (double min, double max) throw()
{
  //  PARALLEL_PRINTF ("line %d min=%f max=%f\n",__LINE__,min,max);
  // simplified variable names

  int mx = nxi_;
  int my = nyi_;
  int m  = mx*my;

  double min2=min;
  double max2=max;
  //  if (! specify_bounds_) {

  min = std::numeric_limits<double>::max();
  max = std::numeric_limits<double>::min();

  // Compute min and max if needed
  //  if (! specify_bounds_) {
  if (image_log_) {
    for (int i=0; i<m; i++) {
      min = MIN(min,log(data_[i]));
      max = MAX(max,log(data_[i]));
    }
  } else {
    for (int i=0; i<m; i++) {
      min = MIN(min,data_[i]);
      max = MAX(max,data_[i]);
    }
  }
  //  }
  //  }
  if (specify_bounds_) {
    min=min2;
    max=max2;
  }
  //  PARALLEL_PRINTF ("line %d min=%f max=%f\n",__LINE__,min,max);

  TRACE1("image_write_() data_ = %p",data_);
  size_t n = map_r_.size();

  // loop over pixels (ix,iy)

  for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = image_log_ ? log(data_[i]) : data_[i];

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
	//	if (value < 0.0) { r=1.0; g=0.0; b=0.0; }

	png_->plot      (ix+1, iy+1, 0.0, 0.0, 0.0);
	png_->plot_blend(ix+1, iy+1, a, r,g,b);

      } else {
	
	// red if out of bounds
	png_->plot(ix+1, iy+1, 1.0, 0.0, 0.0);

      }

      // Plot pixel
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

void OutputImage::reduce_point_ 
(double * data, double value) throw()
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
    *data += value;
    break;
  case reduce_set:
    *data = value;
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_(int ix0, int ix1, int iy0, int iy1, double value)
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
      reduce_point_(&data_[i],value);
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
      reduce_point_(&data_[i],value);
      err += derr;
      if (err >= 0.5) {
	++ix;
	err -= 1.0;
      }
    }
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_x_(int ixm, int ixp, int iy, double value)
{
  ASSERT2("OutputImage::reduce_line_x","! (ixm (%d) <= ixp (%d)",ixm,ixp,ixm<=ixp);
  if (ixp < ixm) { int t = ixp; ixp = ixm; ixm = t; }

  for (int ix=ixm; ix<=ixp; ++ix) {
    int i = ix + nxi_*iy;
    reduce_point_(&data_[i],value);
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_y_(int ix, int iym, int iyp, double value)
{
  ASSERT2("OutputImage::reduce_line_y","! (iym (%d) <= iyp (%d)",iym,iyp,iym<=iyp);
  if (iyp < iym) { int t = iyp; iyp = iym; iym = t; }
  for (int iy=iym; iy<=iyp; ++iy) {
    int i = ix + nxi_*iy;
    reduce_point_(&data_[i],value);
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_box_(int ixm, int ixp, int iym, int iyp, double value)
{
  TRACE5("reduce_box %d %d %d %d %f",ixm,ixp,iym,iyp,value);
  reduce_line_x_(ixm,ixp,iym,value);
  reduce_line_x_(ixm,ixp,iyp,value);
  reduce_line_y_(ixm,iym,iyp,value);
  reduce_line_y_(ixp,iym,iyp,value);
}

//----------------------------------------------------------------------


void OutputImage::reduce_cube_(int ixm, int ixp, int iym, int iyp, double value)
{
  if (ixm > ixp || iym >iyp) return;

  for (int ix=ixm; ix<=ixp; ++ix) {
    for (int iy=iym; iy<=iyp; ++iy) {
      int i = ix + nxi_*iy;
      reduce_point_(&data_[i],value);
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
