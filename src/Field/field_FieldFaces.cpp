// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldFaces.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file field_FieldFaces.cpp
///
/// Detailed description of file field_FieldFaces.cpp

#include "cello.hpp"

#include "field.hpp"

FieldFaces::FieldFaces(FieldBlock * field_block) throw ()
  : array_(),
    index_ghosts_(),
    index_faces_()
{
  allocate_(field_block);
}

//----------------------------------------------------------------------

FieldFaces::~FieldFaces() throw ()
{
  deallocate_();
}

//----------------------------------------------------------------------

FieldFaces::FieldFaces(const FieldFaces & field_faces) throw ()
  : array_(),
    index_ghosts_(),
    index_faces_()
/// @param     FieldFaces  Object being copied
{
  array_        = field_faces.array_;
  index_ghosts_ = field_faces.index_ghosts_;
  index_faces_  = field_faces.index_faces_;

}

//----------------------------------------------------------------------

FieldFaces & FieldFaces::operator= (const FieldFaces & field_faces) throw ()
/// @param     field_faces  Source object of the assignment
///
/// @return    The target assigned object
{
  INCOMPLETE("FieldFaces::operator =");
  return *this;
}

//======================================================================

void FieldFaces::load
(
 const FieldBlock * field_block,
 int field,
 axis_enum axis,
 face_enum face
 ) throw()
{
  const FieldDescr * field_descr = field_block->field_descr();

  if (field == -1) {
    int num_fields = field_descr->field_count();
    // WARNING: recursive
    for (int i=0; i<num_fields; i++) {
      load (field_block,i,axis,face);
    }
  } else {
    // Get precision
    precision_enum precision = field_descr->precision(field);
    // Get field values and face array
    const char * field_values = field_block->field_values(field);
    char * face_values  = &array_[index_faces_[index_(field,axis,face)]];
    // Get field (and face) dimensions
    int n[3];
    field_block->field_size(field,&n[0],&n[1],&n[2]);
    // Compute multipliers for index calculations
    int nd[3] = {1, n[0], n[0]*n[1]};
    // Get ghost depth
    int ng[3];
    field_descr->ghosts(field,&ng[0],&ng[1],&ng[2]);
    // Compute permutation indices
    switch (precision) {
    case precision_single:
      load_precision_((float * )(face_values),
		      (float * )(field_values),
		      n,nd,ng,axis,face);
      break;
    case precision_double:
      load_precision_((double * )(face_values),
		      (double * )(field_values),
		      n,nd,ng,axis,face);
      break;
    case precision_quadruple:
      load_precision_((long double * )(face_values),
		      (long double * )(field_values),
		      n,nd,ng,axis,face);
      break;
    default:
      ERROR("FieldFaces::load", "Unsupported precision");
    }
  }
  
}

//----------------------------------------------------------------------

void FieldFaces::copy
(
 const FieldFaces * field_faces,
 int field,
 axis_enum axis,
 face_enum face
 ) throw()
{
  INCOMPLETE("FieldFaces::copy");
}

//----------------------------------------------------------------------

void FieldFaces::store
(
 FieldBlock * field_block,
 int field,
 axis_enum axis,
 face_enum face
 ) throw()
{
  INCOMPLETE("FieldFaces::store");
}

//----------------------------------------------------------------------

void FieldFaces::allocate_(FieldBlock * field_block) throw()
{
  const FieldDescr * field_descr = field_block->field_descr();

  size_t num_fields = field_descr->field_count();

  index_faces_.reserve (num_fields*3*2);
  index_ghosts_.reserve(num_fields*3*2);

  int array_size = 0;
  int array_index = 0;

  for (size_t index_field = 0; index_field < num_fields; index_field++) {

    // Need element size for alignment adjust below

    precision_enum precision = field_descr->precision(index_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    // Get field block dimensions n3[]
    // Get field_size, which includes ghosts and precision adjustments

    int axis_length[3];
    int field_bytes = field_block->field_size (index_field, 
					      &axis_length[0], 
					      &axis_length[1], 
					      &axis_length[2]);

    // Get ghost depth

    int ghost_count[3];
    field_descr->ghosts(index_field,
			&ghost_count[0],
			&ghost_count[1],
			&ghost_count[2]);

    for (int axis=0; axis<3; axis++) {

      // Get size of one face for given axis

      int face_bytes = ghost_count[axis] * field_bytes / axis_length [axis];

      face_bytes += 
	field_block->adjust_alignment_(face_bytes,bytes_per_element);

      // four arrays of this size: lower and upper, ghosts and faces

      array_size += 4 * face_bytes;

      // Initialize the array_ offsets

      int index_lower = index_(index_field,axis,face_lower);

      index_faces_ [index_lower] = array_index;  array_index += face_bytes;
      index_ghosts_[index_lower] = array_index;  array_index += face_bytes;

      int index_upper = index_(index_field,axis,face_upper);

      index_faces_ [index_upper] = array_index;  array_index += face_bytes;
      index_ghosts_[index_upper] = array_index;  array_index += face_bytes;

    }
  }

  // Allocate the array
  array_.reserve(array_size);

}

//----------------------------------------------------------------------

void FieldFaces::deallocate_() throw()
{
  array_.clear();
  index_ghosts_.clear();
  index_faces_.clear();
}

//----------------------------------------------------------------------

template<class T>
void FieldFaces::load_precision_
(
 T * face_values, T * field_values,
 int n[3], int nd[3], int ng[3],
 axis_enum axis, face_enum face 
)
{
  int iax=(axis+1) % 3;
  int iay=(axis+2) % 3;

  int iz0 = (face == face_lower) ? ng[axis] : n[axis];

  for (int iz = 0; iz <ng[axis]; iz++)  { // 0 <= iz < ng[axis]
    for (int iy=0; iy < n[iay]; iy++) {
      for (int ix=0; ix < n[iax]; ix++) {
	int index_field = ix*nd[iax] + iy*nd[iay] + (iz0+iz)*nd[axis];
	int index_face  = iz + ng[axis]*(ix + n[iax]*iy);
	face_values[index_face] = field_values[index_field];
      }
    }
  }
}

//----------------------------------------------------------------------

template<class T>
void FieldFaces::save_precision_
(
 T * face_values, T * field_values,
 int n[3], int nd[3], int ng[3],
 axis_enum axis, face_enum face 
)
{
  int iax=(axis+1) % 3;
  int iay=(axis+2) % 3;

  int iz0 = (face == face_lower) ? 0 : nd[axis]-ng[axis];

  for (int iz = 0; iz <ng[axis]; iz++)  { // 0 <= iz < ng[axis]
    for (int iy=0; iy < n[iay]; iy++) {
      for (int ix=0; ix < n[iax]; ix++) {
	int index_field = ix*nd[iax] + iy*nd[iay] + (iz0+iz)*nd[axis];
	int index_face  = iz + ng[axis]*(ix + n[iax]*iy);
	field_values[index_field] = face_values[index_face];
      }
    }
  }
}
