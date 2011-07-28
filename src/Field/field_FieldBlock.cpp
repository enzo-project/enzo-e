// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed May 19 18:17:50 PDT 2010
/// @todo     Remove dependence of velocity field names "velocity_[xyz]"
/// @brief    Implementation of the FieldBlock class

#include "cello.hpp"
#include "enzo.hpp"
#include "field.hpp"

#define TEMP_CLEAR_VALUE std::numeric_limits<float>::max() /* in field_FieldBlock.cpp and  mesh_Block.cpp */

//----------------------------------------------------------------------

FieldBlock::FieldBlock ( int nx, int ny, int nz ) throw()
  : array_(0),
    field_values_(),
    ghosts_allocated_(false)
{
  size_[0] = nx;
  size_[1] = ny;
  size_[2] = nz;
}

//----------------------------------------------------------------------

FieldBlock::~FieldBlock() throw()
{  
  deallocate_array();
}

//----------------------------------------------------------------------

FieldBlock::FieldBlock ( const FieldBlock & field_block ) throw ()
{  INCOMPLETE("FieldBlock::FieldBlock"); }

//----------------------------------------------------------------------

FieldBlock & FieldBlock::operator= ( const FieldBlock & field_block ) throw ()
{  INCOMPLETE("FieldBlock::operator=");
  return *this;
}

//----------------------------------------------------------------------

void FieldBlock::size( int * nx, int * ny, int * nz ) const throw()
{
  *nx = size_[0];
  *ny = size_[1];
  *nz = size_[2];
}

//----------------------------------------------------------------------

char * FieldBlock::field_values ( int id_field ) 
  throw (std::out_of_range)
{
  return (unsigned(id_field) < field_values_.size()) ? 
    field_values_.at(id_field) : NULL;
}

//----------------------------------------------------------------------

const char * FieldBlock::field_values ( int id_field ) const 
  throw (std::out_of_range)
{
  return (unsigned(id_field) < field_values_.size()) ? 
    field_values_.at(id_field) : NULL;
}

//----------------------------------------------------------------------

const char * FieldBlock::field_unknowns 
(
 const FieldDescr * field_descr,
 int id_field 
 ) const throw (std::out_of_range)
{
  return (const char *)
    ((FieldBlock *)this) -> field_unknowns(field_descr,id_field);
}

//----------------------------------------------------------------------

char * FieldBlock::field_unknowns 
(
 const FieldDescr * field_descr,
 int id_field 
 ) throw (std::out_of_range)
{
  char * field_unknowns = field_values_.at(id_field);

  if ( ghosts_allocated() ) {

    int gx,gy,gz;
    field_descr->ghosts(id_field,&gx,&gy,&gz);

    bool cx,cy,cz;
    field_descr->centering(id_field,&cx,&cy,&cz);

    int nx,ny,nz;
    size(&nx,&ny,&nz);

    nx += 2*gx + (cx?0:1);
    ny += 2*gy + (cy?0:1);
    nz += 2*gz + (cz?0:1);

    precision_enum precision = field_descr->precision(id_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    field_unknowns += bytes_per_element * (gx + nx*(gy + ny*gz));
  } 

  return field_unknowns;
}

//----------------------------------------------------------------------

void FieldBlock::cell_width(Block * block,
			    double * hx, double * hy, double * hz ) const throw ()
{
  double xm,xp,ym,yp,zm,zp;

  block->lower(&xm,&ym,&zm);
  block->upper(&xp,&yp,&zp);

  if (hx) (*hx) = (xp-xm) / size_[0];
  if (hy) (*hy) = (yp-ym) / size_[1];
  if (hz) (*hz) = (zp-zm) / size_[2];
}

//----------------------------------------------------------------------

// FieldFaces * FieldBlock::field_faces(const FieldDescr * field_descr) throw ()
// {
//   if (field_faces_ == NULL) {
//     field_faces_ = new FieldFaces(field_descr,this);
//   }
//   return field_faces_;
// }

//----------------------------------------------------------------------

void FieldBlock::clear
(
 const FieldDescr * field_descr,
 float              value,
 int                id_field_first,
 int                id_field_last) throw()
{
  if ( array_allocated() ) {
    if (id_field_first == -1) {
      id_field_first = 0;
      id_field_last  = field_descr->field_count() - 1;
    } else if (id_field_last == -1) {
      id_field_last  = id_field_first;
    }
    for (int id_field = id_field_first;
	 id_field <= id_field_last;
	 id_field++) {
      int nx,ny,nz;
      field_size(field_descr,id_field,&nx,&ny,&nz);
      precision_enum precision = field_descr->precision(id_field);
      switch (precision) {
      case precision_single:
	for (int i=0; i<nx*ny*nz; i++) {
	  ((float *)field_values_[id_field])[i] = (float) value;
	}
	break;
      case precision_double:
	for (int i=0; i<nx*ny*nz; i++) {
	  ((double *)field_values_[id_field])[i] = (double) value;
	}
	break;
      case precision_quadruple:
	for (int i=0; i<nx*ny*nz; i++) {
	  ((long double *)field_values_[id_field])[i] = (long double) value;
	}
	break;
      default:
	char buffer[80];
	sprintf (buffer,"Clear called with unsupported precision %s" , 
		 cello::precision_name[precision]);
	ERROR("FieldBlock::clear", buffer);
      }

    }
  } else {
    ERROR("FieldBlock::clear",
		  "Called clear with unallocated arrays");
  }
}

//----------------------------------------------------------------------

void FieldBlock::allocate_array(const FieldDescr * field_descr) throw()
{
  if (! (size_[0] > 0 &&
	 size_[1] > 0 &&
	 size_[2] > 0) ) {
    ERROR ("FieldBlock::allocate_array",
		   "Allocate called with zero field size");
  }

  if ( ! array_allocated() ) {
    
    int padding   = field_descr->padding();
    int alignment = field_descr->alignment();

    int array_size = 0;

    for (int id_field=0; id_field<field_descr->field_count(); id_field++) {

      // Increment array_size, including padding and alignment adjustment

      int nx,ny,nz;       // not needed

      int size = field_size(field_descr,id_field, &nx,&ny,&nz);

      array_size += adjust_padding_   (size,padding);
      array_size += adjust_alignment_ (size,alignment);

    }

    // Adjust for possible initial misalignment

    array_size += alignment - 1;

    // Allocate the array

    array_ = new char [array_size];

    // Initialize field_begin

    char * field_begin = array_ + align_padding_(array_,alignment);

    int field_offset = 0;

    field_values_.reserve(field_descr->field_count());

    for (int id_field=0; id_field<field_descr->field_count(); id_field++) {

      field_values_.push_back(field_begin + field_offset);

      // Increment array_size, including padding and alignment adjustment

      int nx,ny,nz;       // not needed

      int size = field_size(field_descr,id_field,&nx,&ny,&nz);

      field_offset += adjust_padding_  (size,padding);
      field_offset += adjust_alignment_(size,alignment);
    }

    // check if array_size is too big or too small

    if ( ! ( 0 <= (array_size - field_offset)
	     &&   (array_size - field_offset) < alignment)) {
      ERROR ("FieldBlock::allocate_array",
		     "Code error: array size was computed incorrectly");
    };

  } else {
    ERROR ("FieldBlock::allocate_array",
		   "Allocate called with array already allocated");
  }
}

//----------------------------------------------------------------------

void FieldBlock::deallocate_array () throw()
{
  if ( array_allocated() ) {

    delete [] array_;
    array_ = 0;
    field_values_.clear();

  }
  // if (field_faces_ != 0) {
  //   delete field_faces_;
  //   field_faces_ = 0;
  // }

}

//----------------------------------------------------------------------

bool FieldBlock::array_allocated() const throw()
{
  return array_ != 0;
}

//----------------------------------------------------------------------
	
bool FieldBlock::ghosts_allocated() const throw ()
{
  return ghosts_allocated_;
}

//----------------------------------------------------------------------

void FieldBlock::allocate_ghosts(const FieldDescr * field_descr) throw ()
{
  if (! ghosts_allocated() ) {

    std::vector<char *> old_field_values;
    char *              old_array;

    old_array = array_;
    array_ = 0;

    backup_array_ (field_descr,old_field_values);

    ghosts_allocated_ = true;

    allocate_array(field_descr);

    restore_array_ (field_descr, old_field_values);

    delete [] old_array;

  } else {
    WARNING("FieldBlock::allocate_ghosts",
		    "Allocate called with ghosts already allocated");
  }
}

//----------------------------------------------------------------------

void FieldBlock::deallocate_ghosts(const FieldDescr * field_descr) throw ()
{
  if ( ghosts_allocated() ) {

    std::vector<char *> old_field_values;
    char *              old_array;

    old_array = array_;
    array_ = 0;

    backup_array_ (field_descr,old_field_values);

    ghosts_allocated_ = false;

    allocate_array(field_descr);

    restore_array_ (field_descr, old_field_values);

    delete [] old_array;

  } else {
    WARNING("FieldBlock::deallocate_ghosts",
	    "Function called with ghosts not allocated");
  }
}

//----------------------------------------------------------------------
void FieldBlock::refresh_ghosts(const FieldDescr * field_descr) throw()
{
  INCOMPLETE("FieldBlock::refresh_ghosts");
  if ( ! ghosts_allocated() ) {
    WARNING("FieldBlock::refresh_ghosts",
	    "Called with ghosts not allocated: allocating ghosts");
    allocate_ghosts(field_descr);
  }
}

//----------------------------------------------------------------------

void FieldBlock::split
(
 bool split_x, bool split_y, bool split_z,
 FieldBlock ** new_field_blocks ) throw ()
{
}

//----------------------------------------------------------------------

FieldBlock * FieldBlock::merge
(
 bool merge_x, bool merge_y, bool merge_z,
 FieldBlock ** field_blocks ) throw ()
{
  FieldBlock * new_field_block = 0;
  return new_field_block;
}

//----------------------------------------------------------------------
	
void FieldBlock::read (File * file) throw ()
{
  INCOMPLETE("FieldBlock::read");
}

//----------------------------------------------------------------------

void FieldBlock::write (File * file) const throw ()
{
  INCOMPLETE("FieldBlock::write");
}

//----------------------------------------------------------------------

void FieldBlock::set_field_values 
(
 int    id_field, 
 char * field_values) throw()
{
  INCOMPLETE("FieldBlock::set_field_values");
}

//======================================================================


int FieldBlock::adjust_padding_
(
 int size, 
 int padding) const throw ()
{
  return size + padding;
}

//----------------------------------------------------------------------

int FieldBlock::adjust_alignment_
(
 int size, 
 int alignment) const throw ()
{
  return (alignment - (size % alignment)) % alignment;
}

//----------------------------------------------------------------------

int FieldBlock::align_padding_ (char * start, int alignment) const throw()
{ 
  long unsigned start_long = reinterpret_cast<long unsigned>(start);
  return ( alignment - (start_long % alignment) ) % alignment; 
}

//----------------------------------------------------------------------

int FieldBlock::field_size 
(
 const FieldDescr * field_descr,
 int                id_field,
 int              * nx,
 int              * ny,
 int              * nz
 ) const throw()
{

  // Adjust memory usage due to ghosts if needed

  int  gx,gy,gz;
  if ( ghosts_allocated() ) {
    field_descr->ghosts(id_field,&gx,&gy,&gz);
  } else {
    gx = gy = gz = 0;
  }

  // Adjust memory usage due to field centering if needed

  bool cx,cy,cz;
  field_descr->centering(id_field,&cx,&cy,&cz);

  // Compute array size

  *nx = size_[0] + (1-cx) + 2*gx;
  *ny = size_[1] + (1-cy) + 2*gy;
  *nz = size_[2] + (1-cz) + 2*gz;

  // Return array size in bytes

  precision_enum precision = field_descr->precision(id_field);
  int bytes_per_element = cello::sizeof_precision (precision);
  return (*nx) * (*ny) * (*nz) * bytes_per_element;
}

//----------------------------------------------------------------------

void FieldBlock::print (const FieldDescr * field_descr,
			const char * message,
			double lower[3],
			double upper[3]) const throw()
{
#ifndef CELLO_DEBUG
  return;
#endif

  if ( ! array_allocated() ) {
    PARALLEL_PRINTF("%s FieldBlock %p not allocated\n");
  } else {
    int field_count = field_descr->field_count();
    for (int index_field=0; index_field<field_count; index_field++) {

      int nxd,nyd,nzd;
      field_size(field_descr,index_field,&nxd,&nyd,&nzd);
      int gx,gy,gz;
      field_descr->ghosts(index_field,&gx,&gy,&gz);

      int ixm,iym,izm;
      int ixp,iyp,izp;

      // Exclude ghost zones

       ixm = gx;
       iym = gy;
       izm = gz;

       ixp = nxd - gx;
       iyp = nyd - gy;
       izp = nzd - gz;

      // Include ghost zones

      // ixm = 0;
      // iym = 0;
      // izm = 0;

      // ixp = nxd;
      // iyp = nyd;
      // izp = nzd;


      int nx,ny,nz;

      nx = (ixp-ixm);
      ny = (iyp-iym);
      nz = (izp-izm);

      double hx,hy,hz;

      hx = (upper[0]-lower[0])/(nxd-2*gx);
      hy = (upper[1]-lower[1])/(nyd-2*gy);
      hz = (upper[2]-lower[2])/(nzd-2*gz);

      switch (field_descr->precision(index_field)) {
      case precision_single:
	{
	  float * field = (float * ) field_values(index_field);
	  float min = std::numeric_limits<float>::max();
	  float max = std::numeric_limits<float>::min();
	  double sum = 0.0;
	  for (int iz=izm; iz<izp; iz++) {
	    double z = hz*(iz-gz) + lower[axis_z];
	    for (int iy=iym; iy<iyp; iy++) {
	      double y = hy*(iy-gy) + lower[axis_y];
	      for (int ix=ixm; ix<ixp; ix++) {
		double x = hx*(ix-gx) + lower[axis_x];
		int i = ix + nxd*(iy + nyd*iz);
		min = MIN(min,field[i]);
		max = MAX(max,field[i]);
		sum += field[i];
		if (field[i]==TEMP_CLEAR_VALUE) {
		  printf ("TEMP_CLEAR_VALUE match: %s %d (%g %g %g) (%d %d %d)\n",
			  message,index_field,x,y,z,ix,iy,iz);
		}
#ifdef CELLO_DEBUG_VERBOSE
		printf ("%s %d  %f %f %f  %18.14g\n",message,index_field,x,y,z,field[i]);
#endif
	      }
	    }
	  }
	  double avg = sum / (nx*ny*nz);
	  PARALLEL_PRINTF
	    ("%s [%s] %18.14g\n",
	     message ? message : "",
	     field_descr->field_name(index_field).c_str(),
	     avg);

	}
	break;
      case precision_double:
	{
	  double * field = (double * ) field_values(index_field);
	  double min = std::numeric_limits<double>::max();
	  double max = std::numeric_limits<double>::min();
	  double sum = 0.0;
	  for (int iz=izm; iz<izp; iz++) {
	    double z = hz*(iz-gz) + lower[axis_z];
	    for (int iy=iym; iy<iyp; iy++) {
	      double y = hy*(iy-gy) + lower[axis_y];
	      for (int ix=ixm; ix<ixp; ix++) {
		double x = hx*(ix-gx) + lower[axis_x];
		int i = ix + nxd*(iy + nyd*iz);
		min = MIN(min,field[i]);
		max = MAX(max,field[i]);
		sum += field[i];
		sum += field[i];
		if (field[i]==TEMP_CLEAR_VALUE) {
		  printf ("TEMP_CLEAR_VALUE match: %s %d (%g %g %g) (%d %d %d)\n",
			  message,index_field,x,y,z,ix,iy,iz);
		}
#ifdef CELLO_DEBUG_VERBOSE
		printf ("%s %d  %f %f %f  %18.14g\n",message,index_field,x,y,z,field[i]);
#endif
	      }
	    }
	  }
	  double avg = sum / (nx*ny*nz);
	  PARALLEL_PRINTF
	    ("%s [%s] %18.14g\n",
	     message ? message : "",
	     field_descr->field_name(index_field).c_str(),
	     avg);
	}
	break;
      case precision_quadruple:
	{
	  long double * field = 
	    (long double * ) field_values(index_field);
	  long double min = std::numeric_limits<long double>::max();
	  long double max = std::numeric_limits<long double>::min();
	  long double sum = 0.0;
	  for (int iz=izm; iz<izp; iz++) {
	    for (int iy=iym; iy<iyp; iy++) {
	      for (int ix=ixm; ix<ixp; ix++) {
		int i = ix + nxd*(iy + nyd*iz);
		min = MIN(min,field[i]);
		max = MAX(max,field[i]);
		sum += field[i];
	      }
	    }
	  }
	  long double avg = sum / (nx*ny*nz);
	  PARALLEL_PRINTF
	    ("%s FieldBlock[%p,%s] (%d %d %d) [%18.14g %18.14g %18.14g]\n",
	     message ? message : "",this,
	     field_descr->field_name(index_field).c_str(),
	     nx,ny,nz,min,avg,max);
	}
	break;
      default:
	ERROR("FieldBlock::print", "Unsupported precision");
      }
    }
  }
}

//----------------------------------------------------------------------

void FieldBlock::image 
(
 const FieldDescr * field_descr,
 const char * prefix,
 int cycle,
 int ibx, int iby, int ibz
 ) const throw()
{

  int field_count = field_descr->field_count();
  for (int index_field=0; index_field<field_count; index_field++) {

    int nxd,nyd,nzd;
    field_size(field_descr,index_field,&nxd,&nyd,&nzd);

    int nd3[3] = {nxd,nyd,nzd};
    for (int axis=0; axis<3; axis++) {

      EnzoOutputImage * output = new EnzoOutputImage;

      double r[] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0 };
      double g[] = {0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0 };
      double b[] = {0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0 };

      output->image_set_map(7,r,g,b);

      char axis_name[3] = {'x','y','z'};
      char filename[80];
      int ix = (axis+1) % 3;
      int iy = (axis+2) % 3;

      sprintf (filename,"%s-%s-%04d-%c-%d%d%d.png",
	       prefix,field_descr->field_name(index_field).c_str(),
	       cycle,axis_name[axis],ibx,iby,ibz);
  
      switch (field_descr->precision(index_field)) {
      case precision_single:
	output->image
	  (filename,nd3[ix],nd3[iy], (float *) field_values_[index_field],
	   nxd,nyd,nzd, nxd,nyd,nzd, 0.0,0.0,0.0, axis_enum(axis),reduce_avg,-1.0,1.0);
	break;
      case precision_double:
	output->image
	  (filename,nd3[ix],nd3[iy],(double *) field_values_[index_field],
	   nxd,nyd,nzd, nxd,nyd,nzd, 0.0,0.0,0.0, axis_enum(axis),reduce_avg,-1.0,1.0);
	break;
      case precision_quadruple:
	output->image
	  (filename,nd3[ix],nd3[iy],(long double *) field_values_[index_field],
	   nxd,nyd,nzd, nxd,nyd,nzd, 0.0,0.0,0.0, axis_enum(axis),reduce_avg,-1.0,1.0);
	break;
      }

      delete output;

    }
  }
}

//----------------------------------------------------------------------

void FieldBlock::backup_array_ 
( const FieldDescr * field_descr,
  std::vector<char *> & old_field_values )
{
  // save old field_values_

  for (int i=0; i<field_descr->field_count(); i++) {
    old_field_values.push_back(field_values_[i]);
  }
  field_values_.clear();

}

//----------------------------------------------------------------------

void FieldBlock::restore_array_ 
( const FieldDescr * field_descr,
  std::vector<char *> & field_values_from) throw (std::out_of_range)
{
  // copy values
  for (int id_field=0; 
       id_field < field_descr->field_count();
       id_field++) {

    // get "to" field size

    int nx2,ny2,nz2;
    field_size(field_descr, id_field, &nx2,&ny2,&nz2);

    // get "from" field size

    ghosts_allocated_ = ! ghosts_allocated_;

    int nx1,ny1,nz1;
    field_size(field_descr, id_field, &nx1,&ny1,&nz1);

    ghosts_allocated_ = ! ghosts_allocated_;

    // determine offsets to unknowns if ghosts allocated

    int offset1 = (nx1-nx2)/2 + nx1* ( (ny1-ny2)/2 + ny1 * (nz1-nz2)/2 );
    offset1 = MAX (offset1, 0);

    int offset2 = (nx2-nx1)/2 + nx2* ( (ny2-ny1)/2 + ny2 * (nz2-nz1)/2 );
    offset2 = MAX (offset2, 0);

    // determine unknowns size

    int nx = MIN(nx1,nx2);
    int ny = MIN(ny1,ny2);
    int nz = MIN(nz1,nz2);

    // adjust for precision

    precision_enum precision = field_descr->precision(id_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    offset1 *= bytes_per_element;
    offset2 *= bytes_per_element;

    // determine array start

    char * array1 = field_values_from.at(id_field) + offset1;
    char * array2 = field_values_.at(id_field)     + offset2;

    // copy values

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  for (int ip=0; ip<bytes_per_element; ip++) {
	    int i1 = ip + bytes_per_element*(ix + nx1*(iy + ny1*iz));
	    int i2 = ip + bytes_per_element*(ix + nx2*(iy + ny2*iz));
	    array2[i2] = array1[i1];
	  }
	}
      }
    }
  }

  field_values_from.clear();
}
