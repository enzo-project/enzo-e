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

FieldFaces::FieldFaces
(
 const FieldDescr * field_descr,
 FieldBlock * field_block) throw ()
  : array_(),
    index_ghosts_(),
    index_faces_()
{
  allocate_(field_descr,field_block);
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
 const FieldDescr * field_descr,
 const FieldBlock * field_block,
 int                field,
 axis_enum          axis,
 face_enum          face
 ) throw()
{

  if (field == -1) {
    int num_fields = field_descr->field_count();
    // WARNING: recursive
    for (int i=0; i<num_fields; i++) {
      load (field_descr, field_block,i,axis,face);
    }
  } else {
    // Get precision
    precision_enum precision = field_descr->precision(field);
    // Get field values and face array
    const char * field_face = field_block->field_values(field);
    char * array_face  = &array_[index_faces_[index_(field,axis,face)]];
    // Get field (and face) dimensions
    int nd3[3];
    field_block->field_size(field_descr,field,&nd3[0],&nd3[1],&nd3[2]);
    // Compute multipliers for index calculations
    int md3[3] = {1, nd3[0], nd3[0]*nd3[1]};
    // Get ghost depth
    int ng3[3];
    field_descr->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);
    // Compute permutation indices
    switch (precision) {
    case precision_single:
      load_precision_((float * )      (array_face),
		      (const float * )(field_face),
		      nd3,md3,ng3,axis,face);
      break;
    case precision_double:
      load_precision_((double * )      (array_face),
		      (const double * )(field_face),
		      nd3,md3,ng3,axis,face);
      break;
    case precision_quadruple:
      load_precision_((long double * )      (array_face),
		      (const long double * )(field_face),
		      nd3,md3,ng3,axis,face);
      break;
    default:
      ERROR("FieldFaces::load", "Unsupported precision");
    }
  }
  
}

//----------------------------------------------------------------------

void FieldFaces::copy
(
 const FieldDescr * field_descr,
 const FieldFaces * neighbor_faces,
 const FieldBlock * field_block, // required for field_descr and array size
 int field,
 axis_enum axis,
 face_enum face
 ) throw()
{
  if (field == -1) {
    int num_fields = field_descr->field_count();
    // WARNING: recursive
    for (int i=0; i<num_fields; i++) {
      copy (field_descr,neighbor_faces,field_block,i,axis,face);
    }
  } else {
    // Get precision
    precision_enum precision = field_descr->precision(field);
    // Get face array and ghosts array
    int index_face  = index_(field,axis,1-face);
    int index_ghost = index_(field,axis,face);
    const char * array_face = &neighbor_faces->array_[index_faces_[index_face]];
    char *       array_ghost=          &this->array_[index_ghosts_[index_ghost]];
    // Get ghost depth
    int ng3[3];
    field_descr->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);
    // Get field (and face) dimensions
    int nd3[3];
    field_block->field_size(field_descr,field,&nd3[0],&nd3[1],&nd3[2]);
    // Compute permutation indices
    int iax=(axis+1) % 3;
    int iay=(axis+2) % 3;
    // Compute array length
    int n = ng3[axis]*nd3[iax]*nd3[iay];

    switch (precision) {
    case precision_single:
      copy_precision_(      (float * )(array_ghost),
		      (const float * )(array_face),
		      n);
      break;
    case precision_double:
      copy_precision_(      (double * )(array_ghost),
		      (const double * )(array_face),
		      n);
      break;
    case precision_quadruple:
      copy_precision_(      (long double * )(array_ghost),
		      (const long double * )(array_face),
		      n);
      break;
    default:
      ERROR("FieldFaces::store", "Unsupported precision");
    }
  }
}

//----------------------------------------------------------------------

void FieldFaces::store
(
 const FieldDescr * field_descr,
 FieldBlock * field_block,
 int          field,
 axis_enum    axis,
 face_enum    face
 ) throw()
{

  if (field == -1) {
    int num_fields = field_descr->field_count();
    // WARNING: recursive
    for (int i=0; i<num_fields; i++) {
      store (field_descr, field_block,i,axis,face);
    }
  } else {
    // Get precision
    precision_enum precision = field_descr->precision(field);
    // Get field values and face array
    const char * field_ghost = field_block->field_values(field);
    char * array_ghost  = &array_[index_ghosts_[index_(field,axis,face)]];
    // Get field (and face) dimensions
    int nd3[3];
    field_block->field_size(field_descr,field,&nd3[0],&nd3[1],&nd3[2]);
    // Compute multipliers for index calculations
    int md3[3] = {1, nd3[0], nd3[0]*nd3[1]};
    // Get ghost depth
    int ng3[3];
    field_descr->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);
    // Compute permutation indices
    switch (precision) {
    case precision_single:
      store_precision_((float * )      (field_ghost),
		       (const float * )(array_ghost),
		       nd3,md3,ng3,axis,face);
      break;
    case precision_double:
      store_precision_((double * )      (field_ghost),
		       (const double * )(array_ghost),
		       nd3,md3,ng3,axis,face);
      break;
    case precision_quadruple:
      store_precision_((long double * )      (field_ghost),
		       (const long double * )(array_ghost),
		       nd3,md3,ng3,axis,face);
      break;
    default:
      ERROR("FieldFaces::store", "Unsupported precision");
    }
  }
}

//----------------------------------------------------------------------

void FieldFaces::allocate_
(
 const FieldDescr * field_descr,
 FieldBlock * field_block) throw()
{

  size_t num_fields = field_descr->field_count();

  index_faces_.reserve (num_fields*3*2);
  index_ghosts_.reserve(num_fields*3*2);

  int array_size = 0;
  int array_index = 0;

  for (size_t index_field = 0; index_field < num_fields; index_field++) {

    // Need element size for alignment adjust below

    precision_enum precision = field_descr->precision(index_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    // Get field block dimensions nd3[]
    // Get field_size, which includes ghosts and precision adjustments

    int axis_length[3];
    int field_bytes = field_block->field_size 
      (field_descr,index_field, 
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
 T *       array_face, 
 const T * field_face,
 int       nd3[3], 
 int       md3[3], 
 int       ng3[3],
 axis_enum axis, 
 face_enum face 
) throw()
{
  int iax=(axis+1) % 3;
  int iay=(axis+2) % 3;

  int iz0 = (face == face_lower) ? ng3[axis] : nd3[axis]-2*ng3[axis];

  int max=0;
  for (int iz = 0; iz <ng3[axis]; iz++)  { // 0 <= iz < ng3[axis]
    for (int iy=0; iy < nd3[iay]; iy++) {
      for (int ix=0; ix < nd3[iax]; ix++) {
	int index_array_face  = iz + ng3[axis]*(ix + nd3[iax]*iy);
	int index_field_face = ix*md3[iax] + iy*md3[iay] + (iz0+iz)*md3[axis];
	array_face[index_array_face] = 
	  field_face[index_field_face];
	max = MAX(max,index_field_face);
      }
    }
  }
}

//----------------------------------------------------------------------

template<class T>

void FieldFaces::copy_precision_
(
 T *       array_ghost, 
 const T * array_face,
 int       n
) throw()
{
  for (int i = 0; i<n; i++) {
    array_ghost[i] = array_face[i];
  }
}

//----------------------------------------------------------------------

template<class T>
void FieldFaces::store_precision_
(
 T *       field_ghost,
 const T * array_ghost,
 int       nd3[3],
 int       md3[3],
 int       ng3[3],
 axis_enum axis,
 face_enum face 
) throw()
{
  int iax=(axis+1) % 3;
  int iay=(axis+2) % 3;

  int iz0 = (face == face_lower) ? 0 : nd3[axis]-ng3[axis];

  int max=0;
  for (int iz = 0; iz <ng3[axis]; iz++)  { // 0 <= iz < ng3[axis]
    for (int iy=0; iy < nd3[iay]; iy++) {
      for (int ix=0; ix < nd3[iax]; ix++) {
	int index_field_ghost = ix*md3[iax] + iy*md3[iay] + (iz0+iz)*md3[axis];
	int index_array_ghost  = iz + ng3[axis]*(ix + nd3[iax]*iy);
	field_ghost[index_field_ghost] = 
	  array_ghost[index_array_ghost];
	max = MAX(max,index_field_ghost);
      }
    }
  }
}
