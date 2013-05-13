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
  for (int i=0; i<3; i++) {
    ghost_[i] = false;
    face_[i] = 0;
    child_[i] = 0;
  }
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
  for (int i=0; i<3; i++) {
    ghost_[i] = false;
    face_[i] = 0;
    child_[i] = 0;
  }
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
  
  for (int i=0; i<3; i++) {
    ghost_[i] = field_face.ghost_[i];
    face_[i] = field_face.face_[i];
    child_[i] = field_face.child_[i];
  }
  prolong_ =  field_face.prolong_;
  restrict_ =  field_face.restrict_;
}
//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
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
#endif
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

    int ndf3[3];
    field_block_->field_size(field_descr_,field,&ndf3[0],&ndf3[1],&ndf3[2]);

    int ng3[3];
    field_descr_->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);
    
    int im3[3],nf3[3];
    load_loop_limits_ (im3,nf3, ndf3,ng3);

    // if (restrict_) {

    //   // restrict field before loading into array

    //   int im = im3[0] + nd3[0]*(im3[1] + nd3[1]*im3[2]);

    //   int n3_c[3]    = { n3_c[0]/2,  n3_c[1]/2,  n3_c[2]/2 };

    //   restrict_->apply((float *)array,           n3_c, n3_c,
    // 		       &((float *)field_face)[im]), nd3,  n3;
      
    // }
    switch (precision) {
    case precision_single:
      index += load_precision_
	( (float * )       (array_face), 
	  (const float * ) (field_face),
	  ndf3,nf3,im3);
      break;
    case precision_double:
      index += load_precision_
	((double * )       (array_face), 
	 (const double * ) (field_face),
	 ndf3,nf3,im3);
      break;
    case precision_quadruple:
      index += load_precision_ 
	((long double * )      (array_face), 
	 (const long double * )(field_face),
	 ndf3,nf3,im3);
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

    int ndf3[3];
    field_block_->field_size(field_descr_,field,&ndf3[0],&ndf3[1],&ndf3[2]);

    int ng3[3];
    field_descr_->ghosts(field,&ng3[0],&ng3[1],&ng3[2]);

    int im3[3],nf3[3];
    store_loop_limits_ (im3,nf3, ndf3,ng3);

    switch (precision) {
    case precision_single:
      index += store_precision_
	((float * )      (field_ghost), 
	 (const float * )(array_ghost),
	 ndf3,nf3,im3);
      break;
    case precision_double:
      index += store_precision_
	((double * )      (field_ghost),
	 (const double * )(array_ghost),
	 ndf3,nf3,im3);
      break;
    case precision_quadruple:
      index += store_precision_
	((long double * )      (field_ghost),
	 (const long double * )(array_ghost),
	 ndf3,nf3,im3);
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


    int ndf3[3];
    int field_bytes = field_block_->field_size 
      (field_descr_,index_field, &ndf3[0], &ndf3[1], &ndf3[2]);

    int ng3[3];
    field_descr_->ghosts(index_field,&ng3[0],&ng3[1],&ng3[2]);

    int n_old = ndf3[0]*ndf3[1]*ndf3[2];

    if (! ghost_[0]) ndf3[0] -= 2*ng3[0];
    if (! ghost_[1]) ndf3[1] -= 2*ng3[1];
    if (! ghost_[2]) ndf3[2] -= 2*ng3[2];

    int n_new = ndf3[0]*ndf3[1]*ndf3[2];

    field_bytes /= n_old;
    field_bytes *= n_new;

    int face_bytes = field_bytes;
    if (face_[0]) face_bytes = (face_bytes * ng3[0]) / ndf3[0];
    if (face_[1]) face_bytes = (face_bytes * ng3[1]) / ndf3[1];
    if (face_[2]) face_bytes = (face_bytes * ng3[2]) / ndf3[2];

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
 T * array, const T * field_face,
 int ndf3[3], int nf3[3],int im3[3]
) throw()
{
   for (int iz=0; iz <nf3[2]; iz++)  {
      int kz = iz+im3[2];
      for (int iy=0; iy < nf3[1]; iy++) {
	int ky = iy+im3[1];
	for (int ix=0; ix < nf3[0]; ix++) {
	  int kx = ix+im3[0];
	  int index_array = ix +   nf3[0]*(iy +   nf3[1] * iz);
	  int index_field = kx + ndf3[0]*(ky + ndf3[1] * kz);
	  array[index_array] = field_face[index_field];
	}
      }
    }
    return (sizeof(T) * nf3[0] * nf3[1] * nf3[2]);

}

//----------------------------------------------------------------------

template<class T>
size_t FieldFace::store_precision_
(
 T * field_ghost, const T * array,
 int ndf3[3], int nf3[3],int im3[3]
 ) throw()
{
  for (int iz=0; iz <nf3[2]; iz++)  {
    int kz = iz+im3[2];
    for (int iy=0; iy < nf3[1]; iy++) {
      int ky = iy+im3[1];
      for (int ix=0; ix < nf3[0]; ix++) {
	int kx = ix+im3[0];
	int index_array = ix +   nf3[0]*(iy +   nf3[1] * iz);
	int index_field = kx + ndf3[0]*(ky + ndf3[1] * kz);
	field_ghost[index_field] = array[index_array];
      }
    }
  }

  return (sizeof(T) * nf3[0] * nf3[1] * nf3[2]);
}

//----------------------------------------------------------------------

void FieldFace::load_loop_limits_
( int im3[3],int nf3[3], int ndf3[3], int ng3[3])
{
  // NOTES: 4:p12

  for (int axis=0; axis<3; axis++) {

    // starting index im3[axis]

    if (face_[axis] == 0) {

      if (ghost_[axis]) {

	im3[axis] =  0;
	nf3[axis] = ndf3[axis];

      } else {

	im3[axis] = ng3[axis];
	nf3[axis]  = ndf3[axis] - 2*ng3[axis];

      }

      if ( prolong_ ) {
	// adjust for child offset
	nf3[axis]  /= 2;;
	im3[axis] += child_[axis] * nf3[axis];
      }
    }

    if (face_[axis] == -1 || face_[axis] == 1) {

	if (face_[axis] == -1) im3[axis] = ng3[axis];   
	if (face_[axis] == +1) im3[axis] = ndf3[axis] - 2*ng3[axis];

      if (restrict_) { // 2*g ghost depth

	if (face_[axis] == 1) im3[axis] -= ng3[axis];
	nf3[axis] = 2*ng3[axis];

      }	else if (prolong_) { // g/2 ghost depth

	if (face_[axis] == 1) im3[axis] += ng3[axis]/2;
	nf3[axis] = (ng3[axis]+1)/2;

      } else {

	nf3[axis] = ng3[axis];

      }
    }

  }
}

//----------------------------------------------------------------------

void FieldFace::store_loop_limits_
( int im3[3],int nf3[3], int ndf3[3], int ng3[3])
{
  // NOTES: 4:p12

  for (int axis=0; axis<3; axis++) {

    if (face_[axis] == 0) {

      if (ghost_[axis]) {

	im3[axis] =  0;
	nf3[axis] = ndf3[axis];

      } else {

	im3[axis] = ng3[axis];
	nf3[axis]  = ndf3[axis] - 2*ng3[axis];

      }

      if ( restrict_ ) {
	nf3[axis]  /= 2;;
	im3[axis] += child_[axis] * nf3[axis];
      }
    }

    if (face_[axis] == -1 || face_[axis] == 1) {

      if (face_[axis] == -1) im3[axis] = 0;
      if (face_[axis] == +1) im3[axis] = ndf3[axis] - ng3[axis];

      nf3[axis] = ng3[axis];

    }

  }
}

