// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-12
/// @brief    Implementation of the FieldFace class

#include "cello.hpp"
#include "field.hpp"
#include "problem_Prolong.hpp"
#include "problem_Restrict.hpp"

//----------------------------------------------------------------------

FieldFace::FieldFace() throw()
  : field_block_(0),
    field_descr_(0),
    array_(),
    restrict_(0),
    prolong_(0)
{
  ghost_[0] = true;
  ghost_[1] = true;
  ghost_[2] = true;
  face_[0] = 0;
  face_[1] = 0;
  face_[2] = 0;
  child_[0] = 0;
  child_[1] = 0;
  child_[2] = 0;
}

//----------------------------------------------------------------------

FieldFace::FieldFace 
(
 FieldBlock * field_block,
 const FieldDescr * field_descr
 ) throw()
  : field_block_(field_block),
    field_descr_((FieldDescr*)field_descr),
    array_(),
    restrict_(0),
    prolong_(0)
{
  ghost_[0] = true;
  ghost_[1] = true;
  ghost_[2] = true;
  face_[0] = 0;
  face_[1] = 0;
  face_[2] = 0;
  child_[0] = 0;
  child_[1] = 0;
  child_[2] = 0;
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
{
  copy_(field_face);
}

//----------------------------------------------------------------------

FieldFace & FieldFace::operator= (const FieldFace & field_face) throw ()
/// @param     field_face  Source object of the assignment
///
/// @return    The target assigned object
{
  copy_(field_face);
  return *this;
}

//----------------------------------------------------------------------

void FieldFace::copy_(const FieldFace & field_face)
{
  array_        = field_face.array_;
  ghost_[0] = field_face.ghost_[0];
  ghost_[1] = field_face.ghost_[1];
  ghost_[2] = field_face.ghost_[2];
  face_[0] = field_face.face_[0];
  face_[1] = field_face.face_[1];
  face_[2] = field_face.face_[2];
  child_[0] = field_face.child_[0];
  child_[1] = field_face.child_[1];
  child_[2] = field_face.child_[2];
  prolong_ =  field_face.prolong_;
  restrict_ =  field_face.restrict_;
}
//----------------------------------------------------------------------

void FieldFace::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  bool up = p.isUnpacking();

  if (up) field_descr_ = new FieldDescr;
  p | *field_descr_;
  if (up) field_block_ = new FieldBlock;
  p | *field_block_;
  p | array_;
  PUParray(p,face_,3);
  PUParray(p,ghost_,3);
  PUParray(p,child_,3);
  p | *restrict_;  // PUPable
  p | *prolong_;  // PUPable
}

//======================================================================

void FieldFace::load ( int * n, char ** array) throw()
{
  if (array_.size() == 0)  allocate ();

  size_t num_fields = field_descr_->field_count();

  size_t index = 0;

  for (size_t field=0; field< num_fields; field++) {
  
    precision_type precision = field_descr_->precision(field);

    const char * field_face = field_block_->field_values(field);

    char * array_face  = &array_[index];

    int nd3[3];
    field_block_->field_size(field_descr_,field,&nd3[0],&nd3[1],&nd3[2]);

    int ng3[3];
    field_descr_->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);
    
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

    precision_type precision = field_descr_->precision(field);

    char * field_ghost = field_block_->field_values(field);
    
    char * array_ghost  = array + index;

    int nd3[3];
    field_block_->field_size(field_descr_,field,&nd3[0],&nd3[1],&nd3[2]);

    int ng3[3];
    field_descr_->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);

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

  for (size_t index_field = 0; index_field < num_fields; index_field++) {

    // Need element size for alignment adjust below

    precision_type precision = field_descr_->precision(index_field);
    int bytes_per_element = cello::sizeof_precision (precision);


    int nd3[3];
    int field_bytes = field_block_->field_size 
      (field_descr_,index_field, &nd3[0], &nd3[1], &nd3[2]);

    int ng3[3];
    field_descr_->ghosts(index_field,&ng3[0],&ng3[1],&ng3[2]);

    int n_old = nd3[0]*nd3[1]*nd3[2];

    if (! ghost_[0]) nd3[0] -= 2*ng3[0];
    if (! ghost_[1]) nd3[1] -= 2*ng3[1];
    if (! ghost_[2]) nd3[2] -= 2*ng3[2];

    int n_new = nd3[0]*nd3[1]*nd3[2];

    field_bytes /= n_old;
    field_bytes *= n_new;

    int face_bytes = field_bytes;
    if (face_[0]) face_bytes = (face_bytes * ng3[0]) / nd3[0];
    if (face_[1]) face_bytes = (face_bytes * ng3[1]) / nd3[1];
    if (face_[2]) face_bytes = (face_bytes * ng3[2]) / nd3[2];

    face_bytes += 
      field_block_->adjust_alignment_(face_bytes,bytes_per_element);

    array_size += face_bytes;

  }


  array_.resize(array_size);

  for (int i=0; i<array_size; i++) array_[i] = 0;

  return &array_[0];
}

//----------------------------------------------------------------------

void FieldFace::deallocate() throw()
{  array_.clear(); }

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
  bool load;
  int i0[3],n[3];
  loop_limits_ (i0,n, nd3,ng3, load = true);

  for (int iz=0; iz <n[2]; iz++)  {
    int kz = iz+i0[2];
    for (int iy=0; iy < n[1]; iy++) {
      int ky = iy+i0[1];
      for (int ix=0; ix < n[0]; ix++) {
	int kx = ix+i0[0];
	int index_array = ix +   n[0]*(iy +   n[1] * iz);
	int index_field = kx + nd3[0]*(ky + nd3[1] * kz);
	array[index_array] = field_face[index_field];
      }
    }
  }

  return (sizeof(T) * n[0] * n[1] * n[2]);
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
  bool load;
  int i0[3],n[3];
  loop_limits_ (i0,n, nd3,ng3, load = false);

  for (int iz=0; iz <n[2]; iz++)  {
    int kz = iz+i0[2];
    for (int iy=0; iy < n[1]; iy++) {
      int ky = iy+i0[1];
      for (int ix=0; ix < n[0]; ix++) {
	int kx = ix+i0[0];
	int index_array = ix +   n[0]*(iy +   n[1] * iz);
	int index_field = kx + nd3[0]*(ky + nd3[1] * kz);
	field_ghost[index_field] = array[index_array];
      }
    }
  }

  return (sizeof(T) * n[0] * n[1] * n[2]);
}

//----------------------------------------------------------------------

void FieldFace::loop_limits_
( int i0[3],int n[3], int nd3[3], int ng3[3],  bool load )
{
  // NOTES: 4:p12

  bool store = ! load;

  for (int axis=0; axis<3; axis++) {

    // starting index i0[axis]

    if (face_[axis] == 0) {

      if (ghost_[axis]) {

	i0[axis] =  0;
	n[axis] = nd3[axis];

      } else {

	i0[axis] = ng3[axis];
	n[axis]  = nd3[axis] - 2*ng3[axis];

      }

      if ( (load && prolong_) ||
	   (store && restrict_)) {
	// adjust for child offset
	INCOMPLETE("FieldFace::loop_limits_()");
      }
    }

    if (face_[axis] == -1 || face_[axis] == 1) {

      if (load) { // face data
	if (face_[axis] == -1) i0[axis] = ng3[axis];   
	if (face_[axis] == +1) i0[axis] = nd3[axis] - 2*ng3[axis];
      }

      if (store) { // ghost data
	if (face_[axis] == -1) i0[axis] = 0;
	if (face_[axis] == +1) i0[axis] = nd3[axis] - ng3[axis];
      }

      n[axis] = ng3[axis];

      if (load && restrict_) { 
	// adjust for 2g depth
	INCOMPLETE("FieldFace::loop_limits_()");
      }
      if (load && prolong_) {
	// adjust for 1/2 g depth
	INCOMPLETE("FieldFace::loop_limits_()");
      }
    }

  }
}

