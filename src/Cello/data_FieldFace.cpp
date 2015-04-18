// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-12
/// @brief    Implementation of the FieldFace class

#include "cello.hpp"
#include "data.hpp"
#include "problem_Prolong.hpp"
#include "problem_Restrict.hpp"

enum enum_op {
  op_unknown,
  op_load,
  op_store
};

//----------------------------------------------------------------------

FieldFace::FieldFace() throw()
  : field_data_(0),
    array_(),
    restrict_(0),
    prolong_(0),
    field_list_()
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
 FieldData * field_data
 ) throw()
  : field_data_(field_data),
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
  : field_data_(field_face.field_data_),
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
  prolong_  =  field_face.prolong_;
  restrict_ =  field_face.restrict_;
  field_list_ = field_face.field_list_;
}
//----------------------------------------------------------------------

void FieldFace::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  const bool up = p.isUnpacking();

  if (up) field_data_ = new FieldData;
  p | *field_data_;
  p | array_;
  PUParray(p,face_,3);
  PUParray(p,ghost_,3);
  PUParray(p,child_,3);
  p | *restrict_;  // PUPable
  p | *prolong_;  // PUPable
  p | field_list_;
}

//======================================================================

void FieldFace::load ( int * n, char ** array) throw()
{
  if (array_.size() == 0)  allocate ();

  ASSERT("FieldFace::load()",
	 "array_.size() must be > 0",
	 array_.size() > 0);

  const size_t num_fields = field_data_->field_count();

  size_t index_array = 0;

  if (field_list_.size() == 0) {
    field_list_.resize(num_fields);
    for (size_t i=0; i<num_fields; i++) field_list_[i] = i;
  }

  ASSERT("FieldFace::load()",
	 "field_list.size() must be > 0",
	 field_list_.size() > 0);

  for (size_t index_field_list=0;
       index_field_list < field_list_.size();
       index_field_list++) {

    size_t index_field = field_list_[index_field_list];
  
    precision_type precision = field_data_->precision(index_field);

    const void * field_face = field_data_->values(index_field);

    void * array_face  = &array_[index_array];

    int nd3[3],ng3[3],im3[3],n3[3];

    field_data_->field_size(index_field,&nd3[0],&nd3[1],&nd3[2]);
    field_data_->ghosts(index_field,&ng3[0],&ng3[1],&ng3[2]);

    load_loop_limits_ (im3,n3, nd3,ng3);

    if (restrict_) {

      // Restrict field to array

      int nc3[3] = { (n3[0]+1)/2, (n3[1]+1)/2,(n3[2]+1)/2 };

      int im3_array[3] = {0,0,0};

      index_array += restrict_->apply
	(precision, 
	 array_face,nc3,im3_array,nc3, 
	 field_face,nd3,im3,      n3);

    } else {

      // Copy field to array
      switch (precision) {
      case precision_single:
	{
	  float *       array = (float *) array_face;
	  const float * field = (const float *) field_face;
	  index_array += load_ ( array,field, nd3,n3,im3);
	}
	break;
      case precision_double:
	{
	  double *       array = (double *)array_face;
	  const double * field = (const double *) field_face;
	  index_array += load_ ( array,field, nd3,n3,im3);
	}
	break;
      case precision_quadruple:
	{
	  long double *       array = (long double *)array_face;
	  const long double * field = (const long double *) field_face;
	  index_array += load_ ( array,field, nd3,n3,im3);
	}
	break;
      default:
	ERROR("FieldFace::load", "Unsupported precision");
	break;
      }

    }
  }

  *n = array_.size();
  ASSERT("FieldFace::load()",
	 "array size must be > 0",
	 *n > 0);
  *array = &array_[0];

}

//----------------------------------------------------------------------

void FieldFace::store (int n, char * array) throw()
{

  const size_t num_fields = field_data_->field_count();

  size_t index_array = 0;

  for (size_t index_field=0; index_field<num_fields; index_field++) {

    precision_type precision = field_data_->precision(index_field);

    char * field_ghost = field_data_->values(index_field);
    
    char * array_ghost  = array + index_array;

    int nd3[3],ng3[3],im3[3],n3[3];

    field_data_->field_size(index_field,&nd3[0],&nd3[1],&nd3[2]);
    field_data_->ghosts(index_field,&ng3[0],&ng3[1],&ng3[2]);

    store_loop_limits_ (im3,n3, nd3,ng3);

    if (prolong_) {

      // Prolong array to field

      bool need_padding = (ng3[0]%2==1) || (ng3[1]%2==1) || (ng3[2]%2==1);

      ASSERT("FieldFace::store()",
	     "Odd ghost zones not implemented yet: prolong needs padding",
	     ! need_padding);

      int nc3[3] = { (n3[0]+1)/2, (n3[1]+1)/2, (n3[2]+1)/2 };

      int im3_array[3] = {0,0,0};

      index_array += prolong_->apply
	(precision, 
	 field_ghost,nd3,im3,       n3,
	 array_ghost,nc3,im3_array, nc3);

    } else {

      // Copy array to field

      switch (precision) {
      case precision_single:
	{
	  float *       field = (float *)field_ghost;
	  const float * array = (const float *)array_ghost;
	  index_array += store_ (field, array, nd3,n3,im3);
	}
	break;
      case precision_double:
	{
	  double *       field = (double *)field_ghost;
	  const double * array = (const double *)array_ghost;
	  index_array += store_ (field, array, nd3,n3,im3);
	}
	break;
      case precision_quadruple:
	{
	  long double *       field = (long double *)field_ghost;
	  const long double * array = (const long double *)array_ghost;
	  index_array += store_ (field, array, nd3,n3,im3);
	}
	break;
      default:
	ERROR("FieldFace::store", "Unsupported precision");
	break;
      }
    }
  }

  deallocate();
}

//----------------------------------------------------------------------

char * FieldFace::allocate () throw()
{
  size_t num_fields = field_data_->field_count();

  ASSERT("FieldFace::allocate()",
	 "num_fields must be > 0",
	 num_fields > 0);
  
  int array_size = 0;

  for (size_t index_field = 0; index_field < num_fields; index_field++) {

    // Need element size for alignment adjust below

    precision_type precision = field_data_->precision(index_field);
    int bytes_per_element = cello::sizeof_precision (precision);


    int nd3[3];
    int field_bytes = field_data_->field_size 
      (index_field, &nd3[0], &nd3[1], &nd3[2]);

    int ng3[3];
    field_data_->ghosts(index_field,&ng3[0],&ng3[1],&ng3[2]);

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
      field_data_->adjust_alignment_(face_bytes,bytes_per_element);

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
size_t FieldFace::load_
( T * array_face, const T * field_face, 
  int nd3[3], int n3[3],int im3[3] ) throw()
{

  for (int iz=0; iz <n3[2]; iz++)  {
    int kz = iz+im3[2];
    for (int iy=0; iy < n3[1]; iy++) {
      int ky = iy+im3[1];
      for (int ix=0; ix < n3[0]; ix++) {
	int kx = ix+im3[0];
	int index_array = ix +   n3[0]*(iy +   n3[1] * iz);
	int index_field = kx + nd3[0]*(ky + nd3[1] * kz);
	array_face[index_array] = field_face[index_field];
      }
    }
  }

  return (sizeof(T) * n3[0] * n3[1] * n3[2]);

}

//----------------------------------------------------------------------

template<class T> size_t FieldFace::store_
( T * field_ghost, const T * array, int nd3[3], int n3[3],int im3[3] ) throw()
{
  for (int iz=0; iz <n3[2]; iz++)  {
    int kz = iz+im3[2];
    for (int iy=0; iy < n3[1]; iy++) {
      int ky = iy+im3[1];
      for (int ix=0; ix < n3[0]; ix++) {
	int kx = ix+im3[0];
	int index_array = ix +  n3[0]*(iy +  n3[1] * iz);
	int index_field = kx + nd3[0]*(ky + nd3[1] * kz);
	field_ghost[index_field] = array[index_array];
      }
    }
  }

  return (sizeof(T) * n3[0] * n3[1] * n3[2]);
}

//----------------------------------------------------------------------

void FieldFace::load_loop_limits_
( int im3[3],int n3[3], const int nd3[3], const int ng3[3])
{
  // NOTES: 4:p12

  for (int axis=0; axis<3; axis++) {

    // starting index im3[axis]

    if (face_[axis] == 0) {

      if (ghost_[axis]) {

	im3[axis] = 0;
	n3[axis]  = nd3[axis];

      } else {

	im3[axis] = ng3[axis];
	n3[axis]  = nd3[axis] - 2*ng3[axis];

      }

      if ( prolong_ ) {
	// adjust for child offset
	n3[axis] /= 2; 
	im3[axis] += child_[axis] * n3[axis];

	if (ghost_[axis]) {
	  // correct centering for interpolating full coarse block to fine
	  im3[axis] += (1-2*child_[axis]) * ng3[axis]/2 ;
	} else  {
	  //	  n3[axis] += ng3[axis]/2; // should be += ng3[axis]?
	  //	  im3[axis] -= child_[axis]*ng3[axis]/2;
	}

      }

    }

    if (face_[axis] == -1 || face_[axis] == 1) {

      if (face_[axis] == -1) im3[axis] = ng3[axis];   
      if (face_[axis] == +1) im3[axis] = nd3[axis] - 2*ng3[axis];
      n3[axis] = ng3[axis];

      if (restrict_) { // 2*g ghost depth

	if (face_[axis] == 1) im3[axis] -= ng3[axis];
	n3[axis] = 2*ng3[axis];

      }	else if (prolong_) { // g/2 ghost depth

	if (face_[axis] == 1) im3[axis] += ng3[axis]/2;
	n3[axis] = (ng3[axis]+1)/2;

      }
    }
  }

  // Adjust n3 for axes > rank
  if (nd3[0] == 1) n3[0] = 1;
  if (nd3[1] == 1) n3[1] = 1;
  if (nd3[2] == 1) n3[2] = 1;

  check_new_(im3,n3,nd3,ng3,op_load);
}

//----------------------------------------------------------------------

void FieldFace::store_loop_limits_
( int im3[3],int n3[3], const int nd3[3], const int ng3[3])
{
  // NOTES: 4:p12
  for (int axis=0; axis<3; axis++) {

    if (face_[axis] == 0) {

      if (ghost_[axis]) {

	im3[axis] =  0;
	n3[axis] = nd3[axis];

      } else {

	im3[axis] = ng3[axis];
	n3[axis]  = nd3[axis] - 2*ng3[axis];

      }

      if ( restrict_ ) {
	n3[axis] /= 2;
	im3[axis] += child_[axis] * (nd3[axis]-2*ng3[axis])/2;
	if (ghost_[axis]) im3[axis] += ng3[axis]/2;
      }

     
    }

    if (face_[axis] == -1 || face_[axis] == 1) {

      if (face_[axis] == -1) im3[axis] = 0;
      if (face_[axis] == +1) im3[axis] = nd3[axis] - ng3[axis];

      n3[axis] = ng3[axis];

    }

  }
  n3[0] = std::max(n3[0],1);
  n3[1] = std::max(n3[1],1);
  n3[2] = std::max(n3[2],1);

  check_new_(im3,n3,nd3,ng3,op_store);
}

//----------------------------------------------------------------------

void FieldFace::check_new_( int im3[3], int n3[3], const int nd3[3], const int ng3[3],int op_type)
{
  int IM3[3]={im3[0],im3[1],im3[2]};
  int N3[3] ={ n3[0], n3[1], n3[2]};
  int ND3[3]={nd3[0],nd3[1],nd3[2]};
  int NG3[3]={ng3[0],ng3[1],ng3[2]};
  new_loop_limits_ (IM3,N3,ND3,NG3,op_type);

  bool l_im = ((im3[0]==IM3[0])&&(im3[1]==IM3[1])&&(im3[2]==IM3[2]));
  bool l_n  = (( n3[0]==N3[0] )&&( n3[1]==N3[1] )&&( n3[2]==N3[2]));
  bool l_nd = ((nd3[0]==ND3[0])&&(nd3[1]==ND3[1])&&(nd3[2]==ND3[2]));
  bool l_ng = ((ng3[0]==NG3[0])&&(ng3[1]==NG3[1])&&(ng3[2]==NG3[2]));

  static int warning_count = 0;
  const int warning_limit = 100;
  if (warning_count <= warning_limit && ! (l_im && l_n && l_nd && l_ng)) {
    printf ("MISMATCH type %s\n",op_type == op_load ? "load" : "store");
    printf ("MISMATCH im %d %d %d  IM %d %d %d\n",
	    im3[0],im3[1],im3[2],  IM3[0],IM3[1],IM3[2]);
    printf ("MISMATCH  n %d %d %d   N %d %d %d\n",
	     n3[0], n3[1], n3[2],   N3[0], N3[1], N3[2]);
    printf ("MISMATCH nd %d %d %d  ND %d %d %d\n",
	    nd3[0],nd3[1],nd3[2],  ND3[0],ND3[1],ND3[2]);
    printf ("MISMATCH ng %d %d %d  NG %d %d %d\n",
	    ng3[0],ng3[1],ng3[2],  NG3[0],NG3[1],NG3[2]);
    if (prolong_ != NULL) printf ("MISMATCH prolong %p\n",prolong_);
    if (restrict_ != NULL) printf ("MISMATCH restrict %p\n",restrict_);
    printf ("MISMATCH child %d %d %d\n",child_[0],child_[1],child_[2]);
    printf ("MISMATCH ghost %d %d %d\n",ghost_[0],ghost_[1],ghost_[2]);
    printf ("MISMATCH  face %d %d %d\n",face_[0],face_[1],face_[2]);
    printf ("MISMATCH\n");
    ++warning_count;
    if (warning_count == warning_limit) {
      printf ("MISMATCH warning count reached limit of %d\n",
	      warning_count);
    }
  }
}
//----------------------------------------------------------------------

void FieldFace::new_loop_limits_
( int im3[3],int n3[3], const int nd3[3], const int ng3[3], int op_type)
{
  im3[0]=0;
  im3[1]=0;
  im3[2]=0;
  n3[0]=0;
  n3[1]=0;
  n3[2]=0;

  const bool lcopy = ( ! prolong_ && ! restrict_ );

  for (int axis=0; axis<3; axis++) {

    // child offset
    int co = child_[axis]*(nd3[axis]-2*ng3[axis])/2;
    if (lcopy) {
      if (face_[axis] == 0 && ! ghost_[axis]) {
	im3[axis] = ng3[axis];
	n3[axis]  = nd3[axis] - 2*ng3[axis];
      }
      if (face_[axis] == 0 && ghost_[axis]) {
	im3[axis] = 0;
	n3 [axis] = nd3[axis];
      }
      if (face_[axis] == -1 && op_type == op_load) {
	im3[axis] = ng3[axis];
	n3 [axis] = ng3[axis];
      }
      if (face_[axis] == -1 && op_type == op_store) {
	im3[axis] = 0;
	n3 [axis] = ng3[axis];
      }      
      if (face_[axis] == +1 && op_type == op_load) {
	im3[axis] = nd3[axis]-2*ng3[axis];
	n3 [axis] = ng3[axis];
      }
      if (face_[axis] == +1 && op_type == op_store) {
	im3[axis] = nd3[axis]-ng3[axis];
	n3 [axis] = ng3[axis];
      }
    }

    if (prolong_) {
      if (face_[axis] == 0 && ! ghost_[axis] && op_type == op_load) {
	im3[axis] = ng3[axis] + co;
	n3[axis] = (nd3[axis]-2*ng3[axis])/2;
      }
      if (face_[axis] == 0 && ! ghost_[axis] && op_type == op_store) {
	im3[axis] = ng3[axis];
	n3[axis]  = nd3[axis]-2*ng3[axis];
      }	  
      if (face_[axis] == 0 && ghost_[axis] && op_type == op_load) {
	im3[axis] = ng3[axis]/2 + co;
	n3[axis] = nd3[axis]/2;
      }
      if (face_[axis] == 0 && ghost_[axis] && op_type == op_store) {
	im3[axis] = 0;
	n3[axis]  = nd3[axis];
      }
      if (face_[axis] == -1 && op_type == op_load) {
	im3[axis] = ng3[axis];
	n3[axis]  = ng3[axis]/2;
      }
      if (face_[axis] == -1 && op_type == op_store) {
	im3[axis] = 0;
	n3[axis]  = ng3[axis];
      }
      if (face_[axis] == +1 && op_type == op_load) {
	im3[axis] = nd3[axis]-3*ng3[axis]/2;
	n3[axis]  = ng3[axis]/2;
      }
      if (face_[axis] == +1 && op_type == op_store) {
	im3[axis] = nd3[axis]-ng3[axis];
	n3[axis]  = ng3[axis];
      }
    }

    if (restrict_) {
      if (face_[axis] == 0 && !ghost_[axis] && op_type == op_load) {
	im3[axis] = ng3[axis];
	n3[axis]  = nd3[axis]-2*ng3[axis];
      }
      if (face_[axis] == 0 && !ghost_[axis] && op_type == op_store) {
	im3[axis] = ng3[axis] + co;
	n3[axis] = (nd3[axis]-2*ng3[axis])/2;
      }
      if (face_[axis] == 0 && ghost_[axis] && op_type == op_load) {
	im3[axis] = 0;
	n3[axis]  = nd3[axis];
      }
      if (face_[axis] == 0 && ghost_[axis] && op_type == op_store) {
	im3[axis] = ng3[axis]/2 + co;
	n3[axis] = nd3[axis]/2;
      }
      if (face_[axis] == -1 && op_type == op_load) {
	im3[axis] = ng3[axis];
	n3[axis]  = 2*ng3[axis];
      }
      if (face_[axis] == -1 && op_type == op_store) {
	im3[axis] = 0;
	n3[axis]  = ng3[axis];
      }
      if (face_[axis] == +1 && op_type == op_load) {
	im3[axis] = nd3[axis]-3*ng3[axis];
	n3[axis]  = 2*ng3[axis];
      }
      if (face_[axis] == +1 && op_type == op_store) {
	im3[axis] = nd3[axis]-ng3[axis];
	n3[axis]  = ng3[axis];
      }
    }
  }
  n3[0] = std::max(n3[0],1);
  n3[1] = std::max(n3[1],1);
  n3[2] = std::max(n3[2],1);
}
