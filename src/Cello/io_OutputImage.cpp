// See License_CELLO file for license and copyright information

/// @file     io_OutputImage.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar 17 11:14:18 PDT 2011
/// @brief    Implementation of the OutputImage class

#include "cello.hpp"
#include "io.hpp"

// #define TRACE_MEMORY

//----------------------------------------------------------------------
#ifdef TRACE_MEMORY
#  undef TRACE_MEMORY
#  define TRACE_MEMORY(MSG,BYTES) CkPrintf \
  ("TRACE_MEMORY %d : %d %s : %lu\n",                                   \
   CkMyPe(),__LINE__,(MSG),(long unsigned)(BYTES));                     \
  fflush(stdout);
#else
#  define TRACE_MEMORY(MSG,BYTES) /* ... */
#endif

//----------------------------------------------------------------------

OutputImage::OutputImage(int index,
			 const Factory * factory,
			 int process_count,
			 const int root_size_in[3],
			 const int root_blocks_in[3],
			 int min_level, int max_level, int leaf_only,
			 std::string image_type,
			 int image_size[2],
			 std::string image_reduce_type,
			 std::string image_mesh_color,
			 std::string color_particle_attribute,
			 double image_lower[],
			 double image_upper[],
			 int face_rank,
			 int axis,
			 bool image_log,
			 bool image_abs,
			 bool ghost,
                         bool use_min_max,
			 double min_value, double max_value) throw ()
: Output(index,factory),
  image_data_(NULL),
  image_mesh_(NULL),
  color_particle_attribute_(color_particle_attribute),
  axis_(axis),
  use_min_max_(use_min_max),
  min_value_(min_value),
  max_value_(max_value),
  png_(NULL),
  image_type_(image_type),
  face_rank_(face_rank),
  image_log_(image_log),
  image_abs_(image_abs),
  include_ghost_(ghost),
  min_level_(min_level),
  max_level_(max_level),
  leaf_only_(leaf_only)
{
  int root_size[3] =
    {root_size_in[0], root_size_in[1], root_size_in[2]};
  int root_blocks[3] =
    {root_blocks_in[0], root_blocks_in[1], root_blocks_in[2]};
  image_size_[0] = image_size[0];
  image_size_[1] = image_size[1];
  if      (image_reduce_type=="min") { op_reduce_ = reduce_min; }
  else if (image_reduce_type=="max") { op_reduce_ = reduce_max; }
  else if (image_reduce_type=="avg") { op_reduce_ = reduce_avg; }
  else if (image_reduce_type=="sum") { op_reduce_ = reduce_sum; }
  else {
    ERROR1 ("OutputImage::OutputImage()",
	    "Unrecognized output_image_reduce_type %s",
	    image_reduce_type.c_str());
  }

  if      (image_mesh_color=="level")   mesh_color_type_ = mesh_color_level;
  else if (image_mesh_color=="process") mesh_color_type_ = mesh_color_process;
  else if (image_mesh_color=="age")     mesh_color_type_ = mesh_color_age;
  else {
    ERROR1 ("OutputImage::OutputImage()",
	    "Unrecognized output_image_mesh_color %s",
	    image_mesh_color.c_str());
  }

  const int rank = cello::rank();
  if (include_ghost_) {
    for (int i=0; i<rank; i++) {
      root_size[i] *= 2;
    }
  }

  if (image_size_[0] == 0) image_size_[0] = (type_is_mesh_()) ? 2 * root_blocks[0] + 1 : root_size[0];
  if (image_size_[1] == 0) image_size_[1] = (type_is_mesh_()) ? 2* root_blocks[1] + 1 : root_size[1];

  int ngx,ngy,ngz;
  FieldDescr * field_descr = cello::field_descr();
  field_descr->ghost_depth(0,&ngx,&ngy,&ngz);

  if (! type_is_mesh_() && include_ghost_) {
    image_size_[0] += 2*root_blocks[0]*ngx;
    image_size_[1] += 2*root_blocks[1]*ngy;
  }

  // Override default Output::stride_write_: only root writes
  set_stride_write (process_count);
  // Let all processes contribute data when its available
  // (wait stride may be helpful for performance?)
  stride_wait_ = 1;

  // Set default color map to be black and white
  colormap_[0].resize(2);
  colormap_[1].resize(2);
  colormap_[2].resize(2);

  colormap_[0][0] = 0.0;
  colormap_[1][0] = 0.0;
  colormap_[2][0] = 0.0;

  colormap_[0][1] = 1.0;
  colormap_[1][1] = 1.0;
  colormap_[2][1] = 1.0;

  for (int axis=0; axis<3; axis++) {
    image_lower_[axis] = image_lower[axis];
    image_upper_[axis] = image_upper[axis];
  }
}

//----------------------------------------------------------------------

OutputImage::~OutputImage() throw ()
{
  TRACE_MEMORY("delete png_",image_size_[0]*image_size_[1]);
  TRACE_MEMORY("delete image_data_",image_size_[0]*image_size_[1]*sizeof(double));
  TRACE_MEMORY("delete image_mesh_",image_size_[0]*image_size_[1]*sizeof(double));

  delete png_;
  png_ = NULL;
  delete [] image_data_;
  image_data_ = NULL;
  delete [] image_mesh_;
  image_mesh_ = NULL;
}

//----------------------------------------------------------------------

void OutputImage::pup (PUP::er &p)
{
  TRACEPUP;

  Output::pup(p);

  p | colormap_[0];
  p | colormap_[1];
  p | colormap_[2];

  int has_data = (image_data_ != NULL);
  p | has_data;
  if (has_data) {
    if (p.isUnpacking()) {
      TRACE_MEMORY("new image_data",image_size_[0]*image_size_[1]*sizeof(double));
      image_data_ = new double [image_size_[0]*image_size_[1]];
    }
    PUParray(p,image_data_,image_size_[0]*image_size_[1]);
  } else {
    image_data_ = NULL;
  }

  int has_mesh = (image_mesh_ != NULL);
  p | has_mesh;
  if (has_mesh) {
    if (p.isUnpacking()) {
      TRACE_MEMORY("new image_mesh",image_size_[0]*image_size_[1]*sizeof(double));
      image_mesh_ = new double [image_size_[0]*image_size_[1]];
    }
    PUParray(p,image_mesh_,image_size_[0]*image_size_[1]);
  } else {
    image_mesh_ = NULL;
  }
  p | op_reduce_;
  p | mesh_color_type_;
  p | color_particle_attribute_;
  p | axis_;
  p | use_min_max_;
  p | min_value_;
  p | max_value_;
  PUParray(p,image_size_,2);

  WARNING("OutputImage::pup","skipping png");
  // p | *png_;
  if (p.isUnpacking()) png_ = NULL;
  p | image_type_;
  p | face_rank_;
  p | image_log_;
  p | image_abs_;
  p | include_ghost_;
  p | min_level_;
  p | max_level_;
  p | leaf_only_;
  PUParray(p,image_lower_,3);
  PUParray(p,image_upper_,3);
}

//----------------------------------------------------------------------

void OutputImage::set_colormap( std::vector<float>colormap[3])
{
  for (int i=0; i<3; i++) {
    colormap_[i] = colormap[i];
  }
}

//======================================================================

void OutputImage::init () throw()
{
  image_create_();
}

//----------------------------------------------------------------------

void OutputImage::open () throw()
{
  // Open file if writing a single block

  if (is_writer()) {

    std::string file_name = expand_name_ (&file_name_,&file_args_);

    std::string dir_name = directory();

    // Create png object
    Monitor::instance()->print ("Output","writing image file %s",
				(dir_name + "/" + file_name).c_str());
    png_create_(dir_name + "/" + file_name);
    if (chmod (dir_name.c_str(),0755) == -1) {
      ERROR2 ("OutputImage::open()",
	      "chmod() return errno %d: error '%s'",
	      errno,strerror(errno));
    };
  }
}

//----------------------------------------------------------------------

void OutputImage::close () throw()
{
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

void OutputImage::write_block ( const Block *  block ) throw()
// @param block  Block to write
{

  // Exit if Block is not participating in output

  if (! is_active_(block) ) return;

  Field field = ((Data *)block->data())->field();

  const int rank = cello::rank();

  ASSERT("OutputImage::write_block",
	 "axis_ must be = axis_z (2) for 2D problems",
	 ! ( rank == 2 && axis_ != 2));

  ASSERT("OutputImage::write_block",
	 "cannot output rank == 1 problems",
	 ! ( rank == 1) );

  const int IX = (axis_+1) % 3;
  const int IY = (axis_+2) % 3;
  const int IZ = (axis_+0) % 3;

  // Field block size
  int nb3[3];
  field.size(&nb3[0],&nb3[1],&nb3[2]);

  const int level = block->level();

  // Index of (single) field to write

  it_field_index_->first();

  int index_field = (it_field_index_->size() > 0)
    ? it_field_index_->value() : -1;

  // extents of domain
  double dm3[3],dp3[3];
  cello::hierarchy()->lower(dm3,dm3+1,dm3+2);
  cello::hierarchy()->upper(dp3,dp3+1,dp3+2);

  // adjust for image_lower / image_upper

  for (int axis=0; axis<3; axis++) {
    dm3[axis] = std::max(dm3[axis],image_lower_[axis]);
    dp3[axis] = std::min(dp3[axis],image_upper_[axis]);
  }
  // extents of block
  double bm3[3];
  double bp3[3];
  block->lower(bm3,bm3+1,bm3+2);
  block->upper(bp3,bp3+1,bp3+2);

  // image extents of box
  int ixm, ixp;
  int iym, iyp;
  ixm = (bm3[IX]-dm3[IX])/(dp3[IX]-dm3[IX])*image_size_[0];
  iym = (bm3[IY]-dm3[IY])/(dp3[IY]-dm3[IY])*image_size_[1];
  ixp = (bp3[IX]-dm3[IX])/(dp3[IX]-dm3[IX])*image_size_[0];
  iyp = (bp3[IY]-dm3[IY])/(dp3[IY]-dm3[IY])*image_size_[1];

  double h3[3];
  block->cell_width(h3,h3+1,h3+2);

  if (type_is_data_()) {

    if (index_field >= 0) {

      // Get ghost depth

      int ng3[3];
      field.ghost_depth(index_field,&ng3[0],&ng3[1],&ng3[2]);

      // Field array dimensions
      int nd3[3];
      nd3[0] = nb3[0] + 2*ng3[0];
      nd3[1] = nb3[1] + 2*ng3[1];
      nd3[2] = nb3[2] + 2*ng3[2];

      // array offset multipliers
      int d3[3];
      d3[0] = 1;
      d3[1] = nd3[0];
      d3[2] = (rank >= 3) ? nd3[0]*nd3[1] : 0;

      // add block contribution to image

      const char * field_values = (include_ghost_) ?
        field.values(index_field) : field.unknowns(index_field);

      float  * field_float  = (float*)field_values;
      double * field_double = (double*)field_values;

      const int precision = field.precision(index_field);

      double factor = (nb3[IZ] > 1) ? 1.0 / pow(2.0,1.0*level) : 1.0;
      if (rank >= 2 && (std::abs(dm3[IZ] - dp3[IZ]) < h3[IZ])) factor = 1.0;

      int m3[3];
      m3[0] = include_ghost_ ? nd3[0] : nb3[0];
      m3[1] = include_ghost_ ? nd3[1] : nb3[1];
      m3[2] = include_ghost_ ? nd3[2] : nb3[2];
      for (int ix=0; ix<m3[IX]; ix++) {
	double x = bm3[IX] + (ix+0.5)*(bp3[IX]-bm3[IX])/m3[IX];
	int jxm = ixm +  ix   *(ixp-ixm)/m3[IX];
	int jxp = ixm + (ix+1)*(ixp-ixm)/m3[IX]-1;
	if (dm3[IX] <= x && x <= dp3[IX]) {
	  for (int iy=0; iy<m3[IY]; iy++) {
	    double y = bm3[IY] + (iy+0.5)*(bp3[IY]-bm3[IY])/m3[IY];
	    int jym = iym +  iy   *(iyp-iym)/m3[IY];
	    int jyp = iym + (iy+1)*(iyp-iym)/m3[IY]-1;
	    if (dm3[IY] <= y && y <= dp3[IY]) {
	      for (int iz=0; iz<m3[IZ]; iz++) {
		double z = bm3[IZ] + (iz+0.5)*(bp3[IZ]-bm3[IZ])/m3[IZ];
		double zlo = dm3[IZ] - 0.5*h3[IZ];
		double zhi = dp3[IZ] + 0.5*h3[IZ];
		if (rank < 3 || (zlo <= z && z <= zhi)) {
		  int i=ix*d3[IX] + iy*d3[IY] + iz*d3[IZ];
		  double value = 0.0;
		  if (precision == precision_single) {
		    value = field_float[i];
		  } else if (precision == precision_double) {
		    value = field_double[i];
		  }
		  reduce_box_filled_(image_data_,jxm,jxp,jym,jyp, (value*factor));
		}
	      }
	    }
	  }
	}
      }
    }
  }

  if (type_is_mesh_()) {

    // value for mesh
    double value = 0;
    value = mesh_color_(level,block->age());

    if (face_rank_ >= 1) {
      reduce_box_filled_(image_mesh_,ixm,ixp,iym,iyp,value);
    } else {
      reduce_box_filled_(image_mesh_,ixm+3,ixp-3,iym+3,iyp-3,value);
    }
    reduce_box_ (image_mesh_,ixm,ixp,iym,iyp,0.0,reduce_set);

    int xm=(ixm+ixp)/2;
    int ym=(iym+iyp)/2;
    if (face_rank_ <= 1) {
      {
	int if3[3] = {-1,0,0};
	int face_level = block->face_level(if3);
	double face_color = mesh_color_(face_level,0);
	reduce_box_filled_(image_mesh_,ixm+1,ixm+2,ym-1,ym+1, face_color);
      }
      {
	int if3[3] = {1,0,0};
	int face_level = block->face_level(if3);
	double face_color = mesh_color_(face_level,0);
	reduce_box_filled_(image_mesh_,ixp-2,ixp-1,ym-1,ym+1, face_color);
      }
      {
	int if3[3] = {0,-1,0};
	int face_level = block->face_level(if3);
	double face_color = mesh_color_(face_level,0);
	reduce_box_filled_(image_mesh_,xm-1,xm+1,iym+1,iym+2, face_color);
      }
      {
	int if3[3] = {0,1,0};
	int face_level = block->face_level(if3);
	double face_color = mesh_color_(face_level,0);
	reduce_box_filled_(image_mesh_,xm-1,xm+1,iyp-2,iyp-1, face_color);
      }
    }
    if (face_rank_ <= 0) {
      {
	int if3[3] = {-1,-1,0};
	int face_level = block->face_level(if3);
	double face_color = mesh_color_(face_level,0);
	reduce_box_filled_(image_mesh_,ixm+1,ixm+2,iym+1,iym+2, face_color);
      }
      {
	int if3[3] = {1,-1,0};
	int face_level = block->face_level(if3);
	double face_color = mesh_color_(face_level,0);
	reduce_box_filled_(image_mesh_,ixp-2,ixp-1,iym+1,iym+2, face_color);
      }
      {
	int if3[3] = {-1,1,0};
	int face_level = block->face_level(if3);
	double face_color = mesh_color_(face_level,0);
	reduce_box_filled_(image_mesh_,ixm+1,ixm+2,iyp-2,iyp-1, face_color);
      }
      {
	int if3[3] = {1,1,0};
	int face_level = block->face_level(if3);
	double face_color = mesh_color_(face_level,0);
	reduce_box_filled_(image_mesh_,ixp-2,ixp-1,iyp-2,iyp-1, face_color);
      }
    }
  }

  // WRITE PARTICLES

  Particle particle = ((Block * )block)->data()->particle();

  double xdm = dm3[IX];
  double ydm = dm3[IY];
  double xdp = dp3[IX];
  double ydp = dp3[IY];

  for (it_particle_index_->first();
       ! it_particle_index_->done();
       it_particle_index_->next()) {

    const int it = it_particle_index_->value();

    const int ia_color = (color_particle_attribute_ != "") ?
      particle.attribute_index (it, color_particle_attribute_) : -1;

    // get particle attribute stride
    const int ia_x = particle.attribute_position(it,0);
    const int dp = particle.stride(it,ia_x);

    const int nb = particle.num_batches(it);
    const int da = (ia_color != -1) ? particle.stride(it,ia_color) : 0;
    for (int ib=0; ib<nb; ib++) {

      const int np = particle.num_particles(it,ib);
      std::vector<double> position[3];
      position[0].resize(np);
      position[1].resize(np);
      position[2].resize(np);
      particle.position
        (it,ib, position[0].data(),position[1].data(), position[2].data());
      const double * xa = position[IX].data();
      const double * ya = position[IY].data();
      double * pa = (double *) particle.attribute_array (it,ia_color,ib);
      for (int ip=0; ip<np; ip++) {

 	double x = xa[ip*dp];
 	double y = ya[ip*dp];
 	double value = (ia_color == -1) ? 1.0 : pa[ip*da];
	double tx = image_size_[0]*(x - xdm)/(xdp-xdm) - 0.5;
	double ty = image_size_[1]*(y - ydm)/(ydp-ydm) - 0.5;
	int ix0 = floor(tx);
	int iy0 = floor(ty);
	int ix1 = ix0+1;
	int iy1 = iy0+1;
	double ax0 = 1.0 - (tx - floor(tx));
	double ax1 = 1.0 - ax0;
	double ay0 = 1.0 - (ty - floor(ty));
	double ay1 = 1.0 - ay0;

	reduce_point_(image_data_,ix0,iy0,value,ax0*ay0);
	reduce_point_(image_data_,ix1,iy0,value,ax1*ay0);
	reduce_point_(image_data_,ix0,iy1,value,ax0*ay1);
	reduce_point_(image_data_,ix1,iy1,value,ax1*ay1);

      }
    }
  }
}

//----------------------------------------------------------------------

void OutputImage::write_field_data
(
 const FieldData * field_data,
 int index_field) throw()
{
  WARNING("OutputImage::write_field_data",
	  "This function should not be called");
}

//----------------------------------------------------------------------

void OutputImage::write_particle_data
(
 const ParticleData * particle_data,
 int index_particle) throw()
{
  WARNING("OutputImage::write_particle_data",
	  "This function should not be called");
}

//----------------------------------------------------------------------

void OutputImage::prepare_remote (int * n, char ** buffer) throw()
{
  int size = 0;
  int nx = image_size_[0];
  int ny = image_size_[1];

  // Determine buffer size

  size += 2*sizeof(int);        // image_size_[0], image_size_[1]
  size += nx*ny*sizeof(double); // image_data_
  size += nx*ny*sizeof(double); // image_mesh_
  (*n) = size;

  // Allocate buffer (deallocated in cleanup_remote())
  TRACE_MEMORY("new buffer",size);
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

void OutputImage::update_remote  ( int m, char * buffer) throw()
{
  union {
    char   * c;
    double * d;
    int    * i;
  } p ;

  p.c = buffer;

  const int nx = *p.i++;
  const int ny = *p.i++;

  const int n = nx*ny;

  if (op_reduce_ == reduce_min) {
    for (int k=0; k<n; k++) image_data_[k] = std::min(image_data_[k],*p.d++);
    for (int k=0; k<n; k++) image_mesh_[k] = std::min(image_mesh_[k],*p.d++);
  } else if (op_reduce_ == reduce_max) {
    for (int k=0; k<n; k++) image_data_[k] = std::max(image_data_[k],*p.d++);
    for (int k=0; k<n; k++) image_mesh_[k] = std::max(image_mesh_[k],*p.d++);
  } else if (op_reduce_ == reduce_sum) {
    for (int k=0; k<n; k++) image_data_[k] += *p.d++;
    for (int k=0; k<n; k++) image_mesh_[k] += *p.d++;
  } else if (op_reduce_ == reduce_avg) {
    for (int k=0; k<n; k++) image_data_[k] += *p.d++;
    for (int k=0; k<n; k++) image_mesh_[k] += *p.d++;
  } else if (op_reduce_ == reduce_set) {
    for (int k=0; k<n; k++) image_data_[k]  = *p.d++;
    for (int k=0; k<n; k++) image_mesh_[k]  = *p.d++;
  }

}

//----------------------------------------------------------------------

void OutputImage::cleanup_remote  (int * n, char ** buffer) throw()
{
  TRACE_MEMORY("delete buffer",*n);
  delete [] (*buffer);
  (*buffer) = NULL;
}


//======================================================================

double OutputImage::mesh_color_(int level,int age) const
{
  if (mesh_color_type_ == mesh_color_level) {
    return (1.0+level);
  } else if (mesh_color_type_ == mesh_color_process) {
    return (1.0+CkMyPe())/(CkNumPes());
  } else if (mesh_color_type_ == mesh_color_age) {
    return 1.0 / (0.01*age + 1.0);
  } else {
    ERROR1 ("OutputImage::mesh_color_()",
	    "Unknown mesh_color_type_ %d",
	    mesh_color_type_);
    return 0;
  }
}

//----------------------------------------------------------------------

bool OutputImage::is_active_ (const Block * block) const
{
  const int level = block->level();

  if (leaf_only_ && ! block->is_leaf()) {
    return false;
  } else if (! ((min_level_ <= level) && (level <= max_level_))) {
    return false;
  } else {
    return true;
  }
}

//----------------------------------------------------------------------

void OutputImage::png_create_ (std::string filename) throw()
{
  if (is_writer()) {
    const char * file_name = strdup(filename.c_str());
    TRACE_MEMORY("new png_",image_size_[0]*image_size_[1]);
    png_ = new pngwriter(image_size_[0], image_size_[1],0,file_name);
    free ((void *)file_name);
  }
}

//----------------------------------------------------------------------

void OutputImage::png_close_ () throw()
{
  if (is_writer()) {
    png_->close();
    TRACE_MEMORY("delete png_",image_size_[0]*image_size_[1]);
    delete png_;
    png_ = 0;
  }
}

//----------------------------------------------------------------------

void OutputImage::image_create_ () throw()
{
  ASSERT("OutputImage::image_create_",
	 "image_ already created",
	 image_data_ == NULL || image_mesh_ == NULL);

  TRACE_MEMORY("new image_data_",image_size_[0]*image_size_[1]*sizeof(double));
  TRACE_MEMORY("new image_mesh_",image_size_[0]*image_size_[1]*sizeof(double));
  image_data_  = new double [image_size_[0]*image_size_[1]];
  image_mesh_  = new double [image_size_[0]*image_size_[1]];

  const double min = std::numeric_limits<double>::max();
  const double max = -min;

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

  for (int i=0; i<image_size_[0]*image_size_[1]; i++) image_data_[i] = value0;
  for (int i=0; i<image_size_[0]*image_size_[1]; i++) image_mesh_[i] = value0;

}

//----------------------------------------------------------------------

void OutputImage::image_write_ () throw()
{
  // simplified variable names

  int mx = image_size_[0];
  int my = image_size_[1];
  int m  = mx*my;

  double min = std::numeric_limits<double>::max();
  double max = -std::numeric_limits<double>::max();

  // Compute min and max

  if (use_min_max_) {

    min = MIN(min,min_value_);
    max = MAX(max,max_value_);

  } else {

    if (image_log_) {
      for (int i=0; i<m; i++) {
        min = MIN(min,log(fabs(data_(i))));
        max = MAX(max,log(fabs(data_(i))));
      }
    } else if (image_abs_) {
      for (int i=0; i<m; i++) {
        min = MIN(min,fabs(data_(i)));
        max = MAX(max,fabs(data_(i)));
      }
    } else { 
      for (int i=0; i<m; i++) {
        min = MIN(min,data_(i));
        max = MAX(max,data_(i));
      }
    }
  }

  size_t n = colormap_[0].size();

  // loop over pixels (ix,iy)

  for (int ix = 0; ix<mx; ix++) {

    for (int iy = 0; iy<my; iy++) {

      int i = ix + mx*iy;

      double value = data_(i);

      if (image_abs_) value = fabs(value);
      if (image_log_) value = log(fabs(value));

      double r=0.0,g=0.0,b=0.0;

      if (value < min) value = min;
      if (value > max) value = max;


      if (min <= value && value <= max) {

	// map v to lower colormap index
	size_t k =  (n - 1)*(value - min) / (max-min);

	// prevent k == colormap_[0].size()-1, which happens if value == max

	if (k > n - 2) k = n-2;

	// linear interpolate colormap values
	double lo = min +  k   *(max-min)/(n-1);
	double hi = min + (k+1)*(max-min)/(n-1);

	double ratio = (value - lo) / (hi-lo);

	r = (1-ratio)*colormap_[0][k] + ratio*colormap_[0][k+1];
	g = (1-ratio)*colormap_[1][k] + ratio*colormap_[1][k+1];
	b = (1-ratio)*colormap_[2][k] + ratio*colormap_[2][k+1];

	png_->plot      (ix+1, iy+1, r,g,b);

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
  if (type_is_mesh_() && type_is_data_())
    return (image_data_[index] + 0.2*image_mesh_[index])/1.2;
  else if (type_is_data_())
    return image_data_[index];
  else  if (type_is_mesh_())
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

  TRACE_MEMORY("delete image_data",image_size_[0]*image_size_[1]*sizeof(double));
  delete [] image_data_;
  image_data_ = nullptr;

  TRACE_MEMORY("delete image_mesh",image_size_[0]*image_size_[1]*sizeof(double));
  delete [] image_mesh_;
  image_mesh_ = nullptr;
}

//----------------------------------------------------------------------

void OutputImage::reduce_point_
(double * data, int ix, int iy, double value, double alpha) throw()
{
  if ( ! (0 <= ix && ix < image_size_[0])) return;
  if ( ! (0 <= iy && iy < image_size_[1])) return;

  if ( ! (0.0 <= alpha && alpha <= 1.0)) {
    WARNING1 ("OutputImage::reduce_point_()",
	     "Alpha %g is not between 0.0 and 1.0",
	      alpha);
  }
  const int i = ix + image_size_[0]*iy;

  double value_new = 0.0;

  switch (op_reduce_) {
  case reduce_min:
    value_new = alpha*value + (1-alpha)*(data[i]);
    data[i] = std::min(data[i],value_new);
    break;
  case reduce_max:
    value_new = alpha*value + (1-alpha)*(data[i]);
    data[i] = std::max(data[i],value_new);
    break;
  case reduce_avg:
  case reduce_sum:
    value_new = alpha*value;
    data[i] += value_new;
    break;
  case reduce_set:
    value_new = alpha*value + (1-alpha)*(data[i]);
    data[i] = value_new;
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_
(double * data,
 int ix0, int ix1,
 int iy0, int iy1,
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
      reduce_point_(data, ix, iy,value,alpha);
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
      reduce_point_(data,ix,iy,value,alpha);
      err += derr;
      if (err >= 0.5) {
	++ix;
	err -= 1.0;
      }
    }
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_x_
(double * data,
 int ixm, int ixp,
 int iy,
 double value, double alpha)
{
  ASSERT2("OutputImage::reduce_line_x","! (ixm (%d) <= ixp (%d)",ixm,ixp,
          ( ixm <= ixp) );
  if (ixp < ixm) { int t = ixp; ixp = ixm; ixm = t; }

  for (int ix=ixm; ix<=ixp; ++ix) {
    reduce_point_(data,ix,iy,value,alpha);
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_line_y_
(double * data,
 int ix,
 int iym, int iyp,
 double value, double alpha)
{
  ASSERT2("OutputImage::reduce_line_y","! (iym (%d) <= iyp (%d)",iym,iyp,
          ( iym <= iyp ) );
  if (iyp < iym) { int t = iyp; iyp = iym; iym = t; }
  for (int iy=iym; iy<=iyp; ++iy) {
    reduce_point_(data,ix,iy,value,alpha);
  }
}

//----------------------------------------------------------------------

void OutputImage::reduce_box_
(double * data,
 int ixm, int ixp,
 int iym, int iyp,
 double value, reduce_type reduce, double alpha)
{
  reduce_type reduce_save = op_reduce_;
  op_reduce_ = reduce;
  reduce_line_x_(data,ixm,ixp,iym,value,alpha);
  reduce_line_x_(data,ixm,ixp,iyp,value,alpha);
  reduce_line_y_(data,ixm,iym,iyp,value,alpha);
  reduce_line_y_(data,ixp,iym,iyp,value,alpha);
  op_reduce_ = reduce_save;
}

//----------------------------------------------------------------------

void OutputImage::reduce_box_filled_
(double * data,
 int ixm, int ixp,
 int iym, int iyp,
 double value, double alpha)
{
  for (int ix=ixm; ix<=ixp; ++ix) {
    for (int iy=iym; iy<=iyp; ++iy) {
      reduce_point_(data,ix,iy,value,alpha);
    }
  }
}

//----------------------------------------------------------------------
