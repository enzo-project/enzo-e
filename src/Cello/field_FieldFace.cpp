// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-12
/// @brief    Implementation of the FieldFace class

#include "cello.hpp"

#include "field.hpp"

//----------------------------------------------------------------------

FieldFace::FieldFace() throw()
  : field_block_(0),
    field_descr_(0),
    array_()
{
  ghost_[0] = true;
  ghost_[1] = true;
  ghost_[2] = true;
  face_[0] = 0;
  face_[1] = 0;
  face_[2] = 0;
}

//----------------------------------------------------------------------

FieldFace::FieldFace 
(
 FieldBlock * field_block,
 const FieldDescr * field_descr
 ) throw()
  : field_block_(field_block),
    field_descr_((FieldDescr*)field_descr),
    array_()
{
  ghost_[0] = true;
  ghost_[1] = true;
  ghost_[2] = true;
  face_[0] = 0;
  face_[1] = 0;
  face_[2] = 0;

}

//----------------------------------------------------------------------

FieldFace::~FieldFace() throw ()
{
  deallocate();
}

//----------------------------------------------------------------------

FieldFace::FieldFace(const FieldFace & field_face) throw ()
  : field_block_(field_face.field_block_),
    field_descr_(field_face.field_descr_),
    array_()
/// @param     FieldFace  Object being copied
{
  array_        = field_face.array_;
  ghost_[0] = field_face.ghost_[0];
  ghost_[1] = field_face.ghost_[1];
  ghost_[2] = field_face.ghost_[2];
  face_[0] = field_face.face_[0];
  face_[1] = field_face.face_[1];
  face_[2] = field_face.face_[2];
}

//----------------------------------------------------------------------

FieldFace & FieldFace::operator= (const FieldFace & field_face) throw ()
/// @param     field_face  Source object of the assignment
///
/// @return    The target assigned object
{
  array_ = field_face.array_;
  ghost_[0] = field_face.ghost_[0];
  ghost_[1] = field_face.ghost_[1];
  ghost_[2] = field_face.ghost_[2];
  return *this;
}

//======================================================================

void FieldFace::load ( int * n, char ** array) throw()
{
  if (array_.size() == 0) {
    allocate ();
  }

  size_t num_fields = field_descr_->field_count();

  size_t index = 0;

  for (size_t field=0; field< num_fields; field++) {
  
    // Get precision
    precision_type precision = field_descr_->precision(field);

    // Get field values and face array
    const char * field_face = field_block_->field_values(field);

    // Get chunk of array for this field
    char * array_face  = &array_[index];

    // Get field dimensions
    int nd3[3];
    field_block_->field_size(field_descr_,field,&nd3[0],&nd3[1],&nd3[2]);

    // Get ghost depth
    int ng3[3];
    field_descr_->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);

    // Load values for the field
    switch (precision) {
    case precision_single:
      index += load_precision_( (float * )       (array_face), 
				(const float * ) (field_face), nd3,ng3);
      break;
    case precision_double:
      index += load_precision_( (double * )       (array_face), 
				(const double * ) (field_face), nd3,ng3);
      break;
    case precision_quadruple:
      index += load_precision_ ((long double * )      (array_face), 
				(const long double * )(field_face), nd3,ng3);
      break;
    default:
      ERROR("FieldFace::load", "Unsupported precision");
    }
  }

  *n = array_.size();
  *array = &array_[0];
}

//----------------------------------------------------------------------

void FieldFace::store (int n, char * array) throw()
{

  size_t num_fields = field_descr_->field_count();

  size_t index = 0;

  for (size_t field=0; field<num_fields; field++) {

    // Get precision
    precision_type precision = field_descr_->precision(field);

    // Get field values
    char * field_ghost = field_block_->field_values(field);
    
    // Get chunk of array for this field
    char * array_ghost  = array + index;

    // Get field dimensions
    int nd3[3];
    field_block_->field_size(field_descr_,field,&nd3[0],&nd3[1],&nd3[2]);

    // Get ghost depth
    int ng3[3];
    field_descr_->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);

    // Compute permutation indices
    switch (precision) {
    case precision_single:
      index += store_precision_
	((float * )      (field_ghost),
	 (const float * )(array_ghost),
	 nd3,ng3);
      break;
    case precision_double:
      index += store_precision_
	((double * )      (field_ghost),
	 (const double * )(array_ghost),
	 nd3,ng3);
      break;
    case precision_quadruple:
      index += store_precision_
	((long double * )      (field_ghost),
	 (const long double * )(array_ghost),
	 nd3,ng3);
      break;
    default:
      ERROR("FieldFace::store", "Unsupported precision");
    }
  }

  deallocate();
}

//----------------------------------------------------------------------

char * FieldFace::allocate () throw()
{

  size_t num_fields = field_descr_->field_count();

  int array_size = 0;

  // Determine array size

  for (size_t index_field = 0; index_field < num_fields; index_field++) {

    // Need element size for alignment adjust below

    precision_type precision = field_descr_->precision(index_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    // Get field block dimensions nd3[]
    // Get field_size, which includes ghosts and precision adjustments

    int nd3[3];
    int field_bytes = field_block_->field_size 
      (field_descr_,index_field, &nd3[0], &nd3[1], &nd3[2]);

    // Get ghost depth

    int ng3[3];
    field_descr_->ghosts(index_field,&ng3[0],&ng3[1],&ng3[2]);

    // Adjust field_bytes if not storing ghost zones

    int n_old = nd3[0]*nd3[1]*nd3[2];

    if (! ghost_[0]) nd3[0] -= 2*ng3[0];
    if (! ghost_[1]) nd3[1] -= 2*ng3[1];
    if (! ghost_[2]) nd3[2] -= 2*ng3[2];

    int n_new = nd3[0]*nd3[1]*nd3[2];

    field_bytes /= n_old;
    field_bytes *= n_new;

    // Get size of the corner, edge, or face

    int face_bytes = field_bytes;
    if (face_[0]) face_bytes = (face_bytes * ng3[0]) / nd3[0];
    if (face_[1]) face_bytes = (face_bytes * ng3[1]) / nd3[1];
    if (face_[2]) face_bytes = (face_bytes * ng3[2]) / nd3[2];

    face_bytes += 
      field_block_->adjust_alignment_(face_bytes,bytes_per_element);

    // Increment array size

    array_size += face_bytes;

  }

  // Allocate the array

  array_.resize(array_size);
  for (int i=0; i<array_size; i++) array_[i] = 0;

  return &array_[0];
}

//----------------------------------------------------------------------

void FieldFace::deallocate() throw()
{
  array_.clear();
}

//----------------------------------------------------------------------

void FieldFace::loop_limits_
(
 int *ix0,
 int *iy0,
 int *iz0,
 int *nx,
 int *ny,
 int *nz,
 int nd3[3],
 int ng3[3],
 bool load
 )
{
  // offset (ix0,iy0,iz0)

  if (face_[0] == 0) {
    (*ix0) = (ghost_[0]) ? 0 : ng3[0];
  } else {
    if (load) {
      if (face_[0] == -1) (*ix0) = ng3[0];   
      if (face_[0] == +1) (*ix0) = nd3[0] - 2*ng3[0];
    } else {
      if (face_[0] == -1) (*ix0) = 0;
      if (face_[0] == +1) (*ix0) = nd3[0] - ng3[0];
    }
  }
  if (face_[1] == 0) {
    (*iy0) = (ghost_[1]) ? 0 : ng3[1];
  } else {
    if (load) {
      if (face_[1] == -1) (*iy0) = ng3[1];
      if (face_[1] == +1) (*iy0) = nd3[1] - 2*ng3[1];
    } else {
      if (face_[1] == -1) (*iy0) = 0;
      if (face_[1] == +1) (*iy0) = nd3[1] - ng3[1];
    }
  }

  if (face_[2] == 0) {
    (*iz0) = (ghost_[2]) ? 0 : ng3[2];
  } else {
    if (load) {
      if (face_[2] == -1) (*iz0) = ng3[2];
      if (face_[2] == +1) (*iz0) = nd3[2] - 2*ng3[2];
    } else {
      if (face_[2] == -1) (*iz0) = 0;
      if (face_[2] == +1) (*iz0) = nd3[2] - ng3[2];
    }
  }


  if (face_[0] == 0) {
    (*nx) = (ghost_[0]) ? nd3[0] : nd3[0]-2*ng3[0];
  } else {
    (*nx) = ng3[0]; 
  }

  if (face_[1] == 0) {
    (*ny) = (ghost_[1]) ? nd3[1] : nd3[1]-2*ng3[1];
  } else {
    (*ny) = ng3[1]; 
  }

  if (face_[2] == 0) {
    (*nz) = (ghost_[2]) ? nd3[2] : nd3[2]-2*ng3[2];
  } else {
    (*nz) = ng3[2]; 
  }

}

//----------------------------------------------------------------------

template<class T>
size_t FieldFace::load_precision_
(
 T *       array, 
 const T * field_face,
 int       nd3[3], 
 int       ng3[3]
) throw()
{
  int ix0,iy0,iz0;
  int nx,ny,nz;
  bool load;

  // Loop limits

  loop_limits_ (&ix0,&iy0,&iz0,&nx,&ny,&nz, nd3,ng3, load=true);

  TRACE3 (" load face_[0],face_[1],face_[2] = %d %d %d",face_[0],face_[1],face_[2]);
  TRACE3 (" load nd3      = %d %d %d",nd3[0],nd3[1],nd3[2]); 
  TRACE3 (" load ng3      = %d %d %d",ng3[0],ng3[1],ng3[2]); 
  TRACE3 (" load i0       = %d %d %d",ix0,iy0,iz0);
  TRACE3 (" load n        = %d %d %d",nx,ny,nz);
  for (int iz=0; iz <nz; iz++)  {
    int kz = iz+iz0;
    for (int iy=0; iy < ny; iy++) {
      int ky = iy+iy0;
      for (int ix=0; ix < nx; ix++) {
	int kx = ix+ix0;
	int index_array = ix +  nx   *(iy +  ny    * iz);
	int index_field = kx + nd3[0]*(ky + nd3[1] * kz);
	array[index_array] = field_face[index_field];
      }
    }
  }

  return (sizeof(T) * nx * ny * nz);
}

//----------------------------------------------------------------------

template<class T>
size_t FieldFace::store_precision_
(
 T *       field_ghost,
 const T * array,
 int       nd3[3],
 int       ng3[3]
 ) throw()
{

  int ix0,iy0,iz0;
  int nx,ny,nz;
  bool load;

  // Loop limits

  loop_limits_ (&ix0,&iy0,&iz0,&nx,&ny,&nz, nd3,ng3, load=false);

  TRACE3 (" store face_[0],face_[1],face_[2] = %d %d %d",face_[0],face_[1],face_[2]);
  TRACE3 (" store nd3      = %d %d %d",nd3[0],nd3[1],nd3[2]); 
  TRACE3 (" store ng3      = %d %d %d",ng3[0],ng3[1],ng3[2]); 
  TRACE3 (" store i0       = %d %d %d",ix0,iy0,iz0);
  TRACE3 (" store n        = %d %d %d",nx,ny,nz);
  for (int iz=0; iz <nz; iz++)  {
    int kz = iz+iz0;
    for (int iy=0; iy < ny; iy++) {
      int ky = iy+iy0;
      for (int ix=0; ix < nx; ix++) {
	int kx = ix+ix0;
	int index_array = ix +  nx   *(iy +  ny    * iz);
	int index_field = kx + nd3[0]*(ky + nd3[1] * kz);
	field_ghost[index_field] = array[index_array];
      }
    }
  }

  return (sizeof(T) * nx * ny * nz);
}
