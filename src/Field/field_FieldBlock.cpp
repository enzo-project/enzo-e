// $Id: field_FieldBlock.cpp 1388 2010-04-20 23:57:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed May 19 18:17:50 PDT 2010
/// @brief    Implementation of the FieldBlock class

#include "assert.h"

#include "field_FieldBlock.hpp"

FieldBlock::FieldBlock() throw()
  : field_descr_(0),
    array_(0),
    field_values_(),
    ghosts_allocated_(false)
{
  for (int i=0; i<3; i++) {
    dimensions_[i] = 0.0;
    box_lower_[i]  = 0.0;
    box_upper_[i]  = 1.0;
  }
}

//----------------------------------------------------------------------

FieldBlock::~FieldBlock() throw()
{  INCOMPLETE_MESSAGE("FieldBlock::~FieldBlock",""); }

//----------------------------------------------------------------------

FieldBlock::FieldBlock ( const FieldBlock & field_block ) throw ()
{  INCOMPLETE_MESSAGE("FieldBlock::FieldBlock",""); }

//----------------------------------------------------------------------

FieldBlock & FieldBlock::operator= ( const FieldBlock & field_block ) throw ()
{  INCOMPLETE_MESSAGE("FieldBlock::operator=","");
  return *this;
}

//----------------------------------------------------------------------

void FieldBlock::dimensions( int * nx, int * ny, int * nz ) const throw()
{
  *nx = dimensions_[0];
  *ny = dimensions_[1];
  *nz = dimensions_[2];
}

//----------------------------------------------------------------------

void * FieldBlock::field_values ( int id_field ) throw (std::out_of_range)
{
  return field_values_.at(id_field);
}

//----------------------------------------------------------------------

void * FieldBlock::field_unknowns ( int id_field ) throw (std::out_of_range)
{
  char * field_unknowns = field_values_.at(id_field);

  if ( ghosts_allocated() ) {

    int gx,gy,gz;
    field_descr_->ghosts(id_field,&gx,&gy,&gz);

    bool cx,cy,cz;
    field_descr_->centering(id_field,&cx,&cy,&cz);

    int nx,ny,nz;
    dimensions(&nx,&ny,&nz);

    nx += 2*gx + (cx?0:1);
    ny += 2*gy + (cy?0:1);
    nz += 2*gz + (cz?0:1);

    precision_type precision = field_descr_->precision(id_field);
    int bytes_per_element = cello::precision_size (precision);

    field_unknowns += bytes_per_element * (gx + nx*(gy + ny*gz));
  } 

  return field_unknowns;
}

//----------------------------------------------------------------------

void FieldBlock::box_extent
(
 double * lower_x, double * upper_x, 
 double * lower_y, double * upper_y,
 double * lower_z, double * upper_z ) const throw ()
{
  if (lower_x) *lower_x = box_lower_[0];
  if (lower_y) *lower_y = box_lower_[1];
  if (lower_z) *lower_z = box_lower_[2];

  if (upper_x) *upper_x = box_upper_[0];
  if (upper_y) *upper_y = box_upper_[1];
  if (upper_z) *upper_z = box_upper_[2];
}

//----------------------------------------------------------------------

void FieldBlock::cell_width( double * hx, double * hy, double * hz ) const throw ()
{
  *hx = (box_upper_[0] - box_lower_[0]) / dimensions_[0];
  *hy = (box_upper_[1] - box_lower_[1]) / dimensions_[1];
  *hz = (box_upper_[2] - box_lower_[2]) / dimensions_[2];
}

//----------------------------------------------------------------------

FieldDescr * FieldBlock::field_descr() throw ()
{
  return field_descr_;
}

//----------------------------------------------------------------------

void FieldBlock::clear
(
 float value,
 int id_field_first,
 int id_field_last) throw()
{
  if ( array_allocated() ) {
    if (id_field_first == -1) {
      id_field_first = 0;
      id_field_last  = field_descr_->field_count() - 1;
    } else if (id_field_last == -1) {
      id_field_last  = id_field_first;
    }
    for (int id_field = id_field_first;
	 id_field <= id_field_last;
	 id_field++) {
      int nx,ny,nz;
      field_size_(id_field,&nx,&ny,&nz);
      precision_type precision = field_descr_->precision(id_field);
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
	ERROR_MESSAGE("FieldBlock::clear", buffer);
      }

    }
  } else {
    ERROR_MESSAGE("FieldBlock::clear",
		  "Called clear with unallocated arrays");
  }
}

//----------------------------------------------------------------------

void FieldBlock::allocate_array() throw()
{
  if (! (dimensions_[0] > 0 &&
	 dimensions_[1] > 0 &&
	 dimensions_[2] > 0) ) {
    ERROR_MESSAGE ("FieldBlock::allocate_array",
		   "Allocate called with zero field size");
  }

  if ( ! array_allocated() ) {
    
    int padding   = field_descr_->padding();
    int alignment = field_descr_->alignment();

    int array_size = 0;

    for (int id_field=0; id_field<field_descr_->field_count(); id_field++) {

      // Increment array_size, including padding and alignment adjustment

      int nx,ny,nz;       // not needed

      int size = field_size_(id_field, &nx,&ny,&nz);

      array_size += field_size_adjust_(size,padding,alignment);

    }

    // Adjust for possible initial misalignment

    array_size += alignment - 1;

    // Allocate the array

    array_ = new char [array_size];

    // Initialize field_begin

    char * field_begin = array_ + align_padding_(array_,alignment);

    int field_offset = 0;

    for (int id_field=0; id_field<field_descr_->field_count(); id_field++) {

      field_values_.push_back(field_begin + field_offset);

      // Increment array_size, including padding and alignment adjustment

      int nx,ny,nz;       // not needed

      int size = field_size_(id_field,&nx,&ny,&nz);

      field_offset += field_size_adjust_(size,padding,alignment);
    }

    // check if array_size is too big or too small

    if ( ! ( 0 <= (array_size - field_offset)
	     &&   (array_size - field_offset) < alignment)) {
      ERROR_MESSAGE ("FieldBlock::allocate_array",
		     "Code error: array size was computed incorrectly");
    };

  } else {
    ERROR_MESSAGE ("FieldBlock::allocate_array",
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

  } else {
    WARNING_MESSAGE ("FieldBlock::deallocate_array",
		     "Deallocate called with array already deallocated");
  }
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

void FieldBlock::allocate_ghosts() throw ()
{
  if (! ghosts_allocated() ) {

    std::vector<char *> old_field_values;
    char *              old_array;

    old_array = array_;
    array_ = 0;

    backup_array_ (old_field_values);

    ghosts_allocated_ = true;

    allocate_array();

    restore_array_ (old_field_values);

    delete [] old_array;

  } else {
    WARNING_MESSAGE("FieldBlock::allocate_ghosts",
		    "Allocate called with ghosts already allocated");
  }
}

//----------------------------------------------------------------------

void FieldBlock::deallocate_ghosts() throw ()
{
  if ( ghosts_allocated() ) {

    std::vector<char *> old_field_values;
    char *              old_array;

    old_array = array_;
    array_ = 0;

    backup_array_ (old_field_values);

    ghosts_allocated_ = false;

    allocate_array();

    restore_array_ (old_field_values);

    delete [] old_array;

  } else {
    WARNING_MESSAGE("FieldBlock::deallocate_ghosts",
		    "Deallocate called with ghosts already deallocated");
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
	
FieldDescr * FieldBlock::read
(
 File       * file, 
 FieldDescr * field_descr 
 ) throw ()
{
  FieldDescr * new_field_descr = 0;

  return new_field_descr;
}

//----------------------------------------------------------------------

void FieldBlock::write
(
 File       * file,
 FieldDescr * field_descr 
 ) const throw ()
{
  INCOMPLETE_MESSAGE("FieldBlock::write","");
}

//----------------------------------------------------------------------

void FieldBlock::set_dimensions(int nx, int ny, int nz) throw()
{
  if ( ! array_allocated() ) {
    dimensions_[0] = nx;
    dimensions_[1] = ny;
    dimensions_[2] = nz;
  } else {
    // WARNING/ERROR: changing dimensions of allocated array
  }
}

//----------------------------------------------------------------------

void FieldBlock::set_field_values 
(
 int    id_field, 
 char * field_values) throw()
{
  INCOMPLETE_MESSAGE("FieldBlock::set_field_values","");
}

//----------------------------------------------------------------------

void FieldBlock::set_field_descr(FieldDescr * field_descr) throw()
{
  field_descr_ = field_descr;
}

//----------------------------------------------------------------------

void FieldBlock::set_box_extent
(
 double lower_x, double upper_x,
 double lower_y, double upper_y,
 double lower_z, double upper_z ) throw ()

{
  box_lower_[0] = lower_x;
  box_lower_[1] = lower_y;
  box_lower_[2] = lower_z;

  box_upper_[0] = upper_x;
  box_upper_[1] = upper_y;
  box_upper_[2] = upper_z;
}

//======================================================================


int FieldBlock::field_size_adjust_
(
 int size, 
 int padding, 
 int alignment) const throw ()
{
  int field_size = size + padding;
  field_size += (alignment - (field_size % alignment)) % alignment;
  return field_size;
}

//----------------------------------------------------------------------

int FieldBlock::field_size_ 
(
 int id_field,
 int * nx,
 int * ny,
 int * nz
 ) const throw()
{

  // Adjust memory usage due to ghosts if needed

  int  gx,gy,gz;
  if ( ghosts_allocated() ) {
    field_descr_->ghosts(id_field,&gx,&gy,&gz);
  } else {
    gx = gy = gz =0;
  }

  // Adjust memory usage due to field centering if needed

  bool cx,cy,cz;
  field_descr_->centering(id_field,&cx,&cy,&cz);

  // Compute array dimensions

  *nx = dimensions_[0] + (cx ? 0 : 1) + 2*gx;
  *ny = dimensions_[1] + (cy ? 0 : 1) + 2*gy;
  *nz = dimensions_[2] + (cz ? 0 : 1) + 2*gz;

  // Return array size in bytes

  precision_type precision = field_descr_->precision(id_field);
  int bytes_per_element = cello::precision_size (precision);
  return (*nx) * (*ny) * (*nz) * bytes_per_element;
}

//----------------------------------------------------------------------

void FieldBlock::backup_array_ 
( std::vector<char *> & old_field_values )
{
  // save old field_values_

  for (int i=0; i<field_descr_->field_count(); i++) {
    old_field_values.push_back(field_values_[i]);
  }
  field_values_.clear();

}

//----------------------------------------------------------------------

void FieldBlock::restore_array_ 
( std::vector<char *> & field_values_from) throw (std::out_of_range)
{
  // copy values
  for (int id_field=0; 
       id_field < field_descr_->field_count();
       id_field++) {

    // get "to" field size

    int nx2,ny2,nz2;
    field_size_(id_field, &nx2,&ny2,&nz2);

    // get "from" field size

    ghosts_allocated_ = ! ghosts_allocated_;

    int nx1,ny1,nz1;
    field_size_(id_field, &nx1,&ny1,&nz1);

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

    int precision_size = 
      cello::precision_size(field_descr_->precision(id_field));

    offset1 *= precision_size;
    offset2 *= precision_size;

    // determine array start

    char * array1 = field_values_from.at(id_field) + offset1;
    char * array2 = field_values_.at(id_field)     + offset2;

    // copy values

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {
	  for (int ip=0; ip<precision_size; ip++) {
	    int i1 = ip + precision_size*(ix + nx1*(iy + ny1*iz));
	    int i2 = ip + precision_size*(ix + nx2*(iy + ny2*iz));
	    array2[i2] = array1[i1];
	  }
	}
      }
    }
  }

  field_values_from.clear();
}
