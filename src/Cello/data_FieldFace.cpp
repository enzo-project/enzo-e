// See LICENSE_CELLO file for license and copyright information

/// @file     data_FieldFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-04-12
/// @brief    Implementation of the FieldFace class

#include "cello.hpp"
#include "data.hpp"
#include "charm_simulation.hpp"

// #define DEBUG_NEW_REFRESH

long FieldFace::counter[CONFIG_NODE_SIZE] = {0};

#define FORTRAN_NAME(NAME) NAME##_

extern "C" void FORTRAN_NAME(field_face_store_4)
  (float * field, float * array, int * nd3, int * n3);
extern "C" void FORTRAN_NAME(field_face_store_8)
  (double * field, double * array, int * nd3, int * n3);
extern "C" void FORTRAN_NAME(field_face_store_16)
  (long double * field, long double * array, int * nd3, int * n3);

enum enum_op_type {
  op_unknown,
  op_load,
  op_store
};

//----------------------------------------------------------------------

FieldFace::FieldFace 
( const Field & field ) throw()
  :  refresh_type_(refresh_unknown)
{
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d DEBUG FieldFace::FieldFace(Field) %p\n",CkMyPe(),this);
#endif

  ++counter[cello::index_static()];

  for (int i=0; i<3; i++) {
    ghost_[i] = false;
    face_[i]  = 0;
    child_[i] = 0;
  }
}

//----------------------------------------------------------------------

FieldFace::~FieldFace() throw ()
{
  --counter[cello::index_static()];

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d DEBUG FieldFace::~FieldFace(Field) %p\n",CkMyPe(),this);
#endif
}

//----------------------------------------------------------------------

FieldFace::FieldFace(const FieldFace & field_face) throw ()
{
  ++counter[cello::index_static()];

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%d DEBUG FieldFace::FieldFace(FieldFace) %p\n",CkMyPe(),this);
#endif
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

bool FieldFace::operator== (const FieldFace & field_face) throw()
{
  for (int i=0; i<3; i++) {
    if (ghost_[i] != field_face.ghost_[i]) return false;
    if (face_[i]  != field_face.face_[i])  return false;
    if (child_[i] != field_face.child_[i]) return false;
  }
  if (refresh_type_ != field_face.refresh_type_) return false;
  if (field_list_ != field_face.field_list_)     return false;
  return true;
}

//----------------------------------------------------------------------

void FieldFace::copy_(const FieldFace & field_face)
{
  for (int i=0; i<3; i++) {
    ghost_[i] = field_face.ghost_[i];
    face_[i]  = field_face.face_[i];
    child_[i] = field_face.child_[i];
  }
  refresh_type_ = field_face.refresh_type_;
  field_list_ = field_face.field_list_;
  
}
//----------------------------------------------------------------------

void FieldFace::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  PUParray(p,face_,3);
  PUParray(p,ghost_,3);
  PUParray(p,child_,3);
  p | refresh_type_;
  p | field_list_;
  p | restrict_;
  p | prolong_;
}

//======================================================================

void FieldFace::face_to_array ( Field field, int * n, char ** array) throw()
{
  ASSERT("FieldFace::face_to_array()",
	 "field_list.size() must be > 0",
	 field_list_.size() > 0);

  *n = num_bytes_array(field);
  *array = new char [*n];


  ASSERT("FieldFace::face_to_array()",
	 "array size must be > 0",
	 *n > 0);

  face_to_array (field, *array);

}

//----------------------------------------------------------------------
void FieldFace::face_to_array ( Field field,char * array) throw()
{

  size_t index_array = 0;

#ifdef DEBUG_NEW_REFRESH
  print("face_to_array");
#endif
  for (size_t i_f=0; i_f < field_list_.size(); i_f++) {

    size_t index_field = field_list_[i_f];
  
    precision_type precision = field.precision(index_field);

    const void * field_face = field.values(index_field);

    char * array_face  = &array[index_array];

    int nd3[3],ng3[3],im3[3],n3[3];

    field.field_size(index_field,&nd3[0],&nd3[1],&nd3[2]);
    field.ghost_depth(index_field,&ng3[0],&ng3[1],&ng3[2]);

    loop_limits (im3,n3,nd3,ng3,op_load);

    if (refresh_type_ == refresh_coarse) {

      // Restrict field to array

      int nc3[3] = { (n3[0]+1)/2, (n3[1]+1)/2,(n3[2]+1)/2 };

      int im3_array[3] = {0,0,0};

      Simulation * simulation = proxy_simulation.ckLocalBranch();
      Problem * problem   = simulation->problem();
      Restrict * restrict = restrict_ ? restrict_ : problem->restrict();

      index_array += restrict->apply
	(precision, 
	 array_face,nc3,im3_array,nc3, 
	 field_face,nd3,im3,      n3);

    } else {

      // Copy field to array
      switch (precision) {
      case precision_single:
	{
	  index_array += load_ ( (float*)array_face,
				 (const float *) field_face, nd3,n3,im3);
	}
	break;
      case precision_double:
	{
	  index_array += load_ ( (double *)array_face,
				 (const double *) field_face, nd3,n3,im3);
	}
	break;
      case precision_quadruple:
	{
	  index_array += load_ ( (long double *)array_face,
				 (const long double *) field_face, nd3,n3,im3);
	}
	break;
      default:
	ERROR("FieldFace::face_to_array", "Unsupported precision");
	break;
      }

    }
  }

}

//----------------------------------------------------------------------

void FieldFace::array_to_face (char * array, Field field) throw()
{
#ifdef DEBUG_NEW_REFRESH
  print("array_to_face");
#endif
  size_t index_array = 0;

  for (size_t i_f=0; i_f < field_list_.size(); i_f++) {

    size_t index_field = field_list_[i_f];

    precision_type precision = field.precision(index_field);

    char * field_ghost = field.values(index_field);
    
    char * array_ghost  = array + index_array;

    int nd3[3],ng3[3],im3[3],n3[3];

    field.field_size(index_field,&nd3[0],&nd3[1],&nd3[2]);
    field.ghost_depth(index_field,&ng3[0],&ng3[1],&ng3[2]);

    loop_limits (im3,n3,nd3,ng3,op_store);

    if (refresh_type_ == refresh_fine) {

      // Prolong array to field

      bool need_padding = (ng3[0]%2==1) || (ng3[1]%2==1) || (ng3[2]%2==1);

      ASSERT("FieldFace::array_to_face()",
	     "Odd ghost zones not implemented yet: prolong needs padding",
	     ! need_padding);

      int nc3[3] = { (n3[0]+1)/2, (n3[1]+1)/2, (n3[2]+1)/2 };

      int im3_array[3] = {0,0,0};

      Simulation * simulation = proxy_simulation.ckLocalBranch();
      Problem * problem   = simulation->problem();

      Prolong * prolong = prolong_ ? prolong_ : problem->prolong();

      index_array += prolong->apply
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
	ERROR("FieldFace::array_to_face()", "Unsupported precision");
	break;
      }
    }
  }
}

//----------------------------------------------------------------------

void FieldFace::face_to_face (Field field_src, Field field_dst)
{
#ifdef DEBUG_NEW_REFRESH
  CkPrintf("%p face_to_face\n",this);
#endif
  
  int M3[3],G3[3],IS3[3],ID3[3],NS3[3],ND3[3];
  int m3[3],g3[3],is3[3],id3[3],ns3[3],nd3[3];
  
  field_src.field_size (0,&M3[0],&M3[1],&M3[2]);
  field_src.ghost_depth(0,&G3[0],&G3[1],&G3[2]);

  //    invert_face();
  loop_limits (IS3,NS3,M3,G3,op_load);
  invert_face();

  loop_limits (ID3,ND3,M3,G3,op_store);
  invert_face();

  for (size_t i_f=0; i_f < field_list_.size(); i_f++) {

    size_t index_field = field_list_[i_f];

    if (! field_src.is_centered(index_field)) {
      // recompute loop limits if Field is not centered
      field_src.field_size (0,&m3[0],&m3[1],&m3[2]);
      field_src.ghost_depth(0,&g3[0],&g3[1],&g3[2]);
      loop_limits (is3,ns3,m3,g3,op_load);
      invert_face();
      loop_limits (id3,nd3,m3,g3,op_store);
      invert_face();
    }

    precision_type precision = field_src.precision(index_field);

    char * values_src = field_src.values(index_field);
    char * values_dst = field_dst.values(index_field);
    
    Simulation * simulation = proxy_simulation.ckLocalBranch();
    Problem * problem   = simulation->problem();

    if (refresh_type_ == refresh_fine) {

      // Prolong field

      bool need_padding = (g3[0]%2==1) || (g3[1]%2==1) || (g3[2]%2==1);

      ASSERT("FieldFace::array_to_face()",
	     "Odd ghost zones not implemented yet: prolong needs padding",
	     ! need_padding);

      Prolong * prolong = prolong_ ? prolong_ : problem->prolong();

      prolong->apply (precision, 
		      values_dst,m3,id3, nd3,
		      values_src,m3,is3, ns3);

    } else if (refresh_type_ == refresh_coarse) {

      // Restrict field

      Restrict * restrict = restrict_ ? restrict_ : problem->restrict();

      restrict->apply (precision, 
		       values_dst,m3,id3, nd3,
		       values_src,m3,is3, ns3);

    } else {

      // Copy faces to ghosts

      switch (precision) {
      case precision_single:
	{
	  const float * vs = (const float *) values_src;
	  float       * vd = (float *) values_dst;
	  copy_ (vd,m3,nd3,id3, vs,m3,ns3,is3);
	}
	break;
      case precision_double:
	{
	  const double * vs = (const double *) values_src;
	  double       * vd = (double *) values_dst;
	  copy_ (vd,m3,nd3,id3, vs,m3,ns3,is3);
	}
	break;
      case precision_quadruple:
	{
	  const long double * vs = (const long double *) values_src;
	  long double       * vd = (long double *) values_dst;
	  copy_ (vd,m3,nd3,id3, vs,m3,ns3,is3);
	}
	break;
      default:
	ERROR("FieldFace::array_to_face()", "Unsupported precision");
	break;
      }
    }
    if (! field_src.is_centered(index_field) ) {
      // reinitialize loop limits to centered
      m3[0] =  M3[0];  m3[1] =  M3[1];  m3[2] =  M3[2];
      g3[0] =  G3[0];  g3[1] =  G3[1];  g3[2] =  G3[2];
      is3[0] = IS3[0]; is3[1] = IS3[1]; is3[2] = IS3[2];
      id3[0] = ID3[0]; id3[1] = ID3[1]; id3[2] = ID3[2];
      ns3[0] = NS3[0]; ns3[1] = NS3[1]; ns3[2] = NS3[2];
      nd3[0] = ND3[0]; nd3[1] = ND3[1]; nd3[2] = ND3[2];
    }


  }
}

//----------------------------------------------------------------------

int FieldFace::num_bytes_array(Field field) throw()
{
  int array_size = 0;

  for (size_t i_f=0; i_f < field_list_.size(); i_f++) {

    size_t index_field = field_list_[i_f];

    precision_type precision = field.precision(index_field);
    int bytes_per_element = cello::sizeof_precision (precision);

    int nd3[3],ng3[3],im3[3],n3[3];

    field.field_size (index_field,&nd3[0],&nd3[1],&nd3[2]);
    field.ghost_depth(index_field,&ng3[0],&ng3[1],&ng3[2]);

    int op_type = (refresh_type_ == refresh_fine) ? op_load : op_store;
    loop_limits (im3,n3,nd3,ng3,op_type);

    array_size += n3[0]*n3[1]*n3[2]*bytes_per_element;

  }

  ASSERT("FieldFace::num_bytes_array()",
	 "array_size must be > 0, maybe field_list_.size() is 0?",
	 array_size);

  return array_size;

}

//----------------------------------------------------------------------

int FieldFace::data_size () const
{
  int count = 0;

  count += 3*sizeof(int);  // face_[3]
  count += 3*sizeof(bool); // ghost_[3]
  count += 3*sizeof(int);  // child_[3];
  count += 1*sizeof(int);  // refresh_type_ (restrict,prolong,copy)
  count += (1+field_list_.size()) * sizeof(int);

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%s:%d data_size %d\n",__FILE__,__LINE__,count);
  fflush(stdout);
#endif

  return count;

}

//----------------------------------------------------------------------

char * FieldFace::save_data (char * buffer) const
{
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%s:%d save_data %p\n",__FILE__,__LINE__,buffer);
  fflush(stdout);
#endif

  char * p = buffer;
  int n;

  memcpy(p,face_, n=3*sizeof(int));  
  p+=n;

  memcpy(p,ghost_,n=3*sizeof(bool)); 
  p+=n;

  memcpy(p,child_,n=3*sizeof(int));  
  p+=n;

  memcpy(p,&refresh_type_,n=sizeof(int));  
  p+=n;

  int length = field_list_.size();

  memcpy(p,&length,n=sizeof(int)); 
  p+=n;

  memcpy(p,&field_list_[0],n=length*sizeof(int)); 
  p+=n;

  return p;
}

//----------------------------------------------------------------------

char * FieldFace::load_data (char * buffer)
{
#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("load_data\n"); fflush(stdout);
#endif

#ifdef DEBUG_NEW_REFRESH
  CkPrintf ("%s:%d load_data buffer %p\n",__FILE__,__LINE__,buffer);
  fflush(stdout);
#endif

  char * p = buffer;
  int n;

#ifdef DEBUG_NEW_REFRESH
  int *pi = (int *)p;
  CkPrintf ("p = %p\n",p);
  CkPrintf ("p[0] = %d\n",pi[0]);
#endif
  memcpy(face_,p, n=3*sizeof(int));
  p+=n;

  memcpy(ghost_,p,n=3*sizeof(bool));
  p+=n;

  memcpy(child_,p,n=3*sizeof(int));
  p+=n;

  memcpy(&refresh_type_,p,n=sizeof(int));
  p+=n;

  int length;

  memcpy(&length,p,n=sizeof(int)); 
  p+=n;

  field_list_.resize(length);

  memcpy(&field_list_[0],p,n=(field_list_.size()*sizeof(int)));
  p+=n;
  
  return p;
}

//======================================================================

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
( T * ghost, const T * array, int nd3[3], int n3[3],int im3[3] ) throw()
{

#define FORTRAN_STORE

#ifdef FORTRAN_STORE

  // This is to get around a bug on SDSC Comet where this function
  // crashes with -O3 (See Enzo-P / Cello bug report #90)
  // http://client64-249.sdsc.edu/cello-bug/show_bug.cgi?id=90

  union {
    float *       ghost_4;
    double *      ghost_8;
    long double * ghost_16;
  };
  union {
    float *       array_4;
    double *      array_8;
    long double * array_16;
  };

  ghost_4 = (float *) ghost;
  array_4 = (float *) array;
  
  int im = im3[0] + nd3[0]*(im3[1] + nd3[1]*im3[2]);

  if (sizeof(T)==sizeof(float)) {
    FORTRAN_NAME(field_face_store_4)(ghost_4 + im,   array_4, nd3,n3);
  } else if (sizeof(T)==sizeof(double)) {
    FORTRAN_NAME(field_face_store_8)(ghost_8 + im,   array_8, nd3,n3);
  } else if (sizeof(T)==sizeof(long double)) {
    FORTRAN_NAME(field_face_store_16)(ghost_16 + im, array_16, nd3,n3);
  } else {
    ERROR1 ("FieldFace::store_()",
	   "unknown float precision sizeof(T) = %d\n",sizeof(T));
  }

#else

   for (int iz=0; iz <n3[2]; iz++)  {
     int kz = iz+im3[2];
     for (int iy=0; iy < n3[1]; iy++) {
       int ky = iy+im3[1];
       for (int ix=0; ix < n3[0]; ix++) {
   	int kx = ix+im3[0];
   	int index_array = ix +  n3[0]*(iy +  n3[1] * iz);
   	int index_field = kx + nd3[0]*(ky + nd3[1] * kz);
   	ghost[index_field] = array[index_array];
       }
     }
   }

#endif

  return (sizeof(T) * n3[0] * n3[1] * n3[2]);

}

//----------------------------------------------------------------------

template<class T> void FieldFace::copy_
( T       * vd, int md3[3],int nd3[3],int id3[3],
  const T * vs, int ms3[3],int ns3[3],int is3[3]) throw()
{
  for (int iz=0; iz <ns3[2]; iz++)  {
    for (int iy=0; iy < ns3[1]; iy++) {
      for (int ix=0; ix < ns3[0]; ix++) {
	int i_src = (ix+is3[0]) + ms3[0]*((iy+is3[1]) + ms3[1] * (iz+is3[2]));
	int i_dst = (ix+id3[0]) + md3[0]*((iy+id3[1]) + md3[1] * (iz+id3[2]));
      	vd[i_dst] = vs[i_src];
      }
    }
  }
}

//----------------------------------------------------------------------

void FieldFace::loop_limits
( int im3[3],int n3[3], const int nd3[3], const int ng3[3], int op_type)
{
  im3[0]=0;
  im3[1]=0;
  im3[2]=0;
  n3[0]=0;
  n3[1]=0;
  n3[2]=0;

  const bool lcopy = (refresh_type_ == refresh_same);

  for (int axis=0; axis<3; axis++) {

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

    // adjust limits to include ghost zones for oblique edges/corners
    // at coarse-fine level interfaces
    
    const bool full_block = (face_[0] == 0 && face_[1] == 0 && face_[2] == 0);

    // child offset: 0 or n/2

    const int co = child_[axis]*(nd3[axis]-2*ng3[axis])/2;

    if (refresh_type_ == refresh_fine) {

      if (face_[axis] == 0 && ! ghost_[axis] && op_type == op_load) {
	im3[axis] = ng3[axis] + co;
	n3[axis] = (nd3[axis]-2*ng3[axis])/2;

	// Bug #70 fix: always include ghosts in finer block when
	// face_[axis] = 0 see notes 150811

	if (! full_block) {
	  if (child_[axis] == 1) {
	    im3[axis] -= ng3[axis]/2;
	  }
	  n3[axis] += ng3[axis]/2;
	}

      }
      if (face_[axis] == 0 && ghost_[axis] && op_type == op_load) {
	im3[axis] = ng3[axis]/2 + co;
	n3[axis] = nd3[axis]/2;
      }
      if (face_[axis] == 0 && ! ghost_[axis] && op_type == op_store) {
	im3[axis] = ng3[axis];
	n3[axis]  = nd3[axis]-2*ng3[axis];

	// Bug #70 fix: always include ghosts in finer block when
	// face_[axis] = 0 see notes 150811

	if (! full_block) {
	  if (child_[axis] == 1) {
	    im3[axis] -= ng3[axis];
	  }
	  n3[axis] += ng3[axis];
	}

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

    if (refresh_type_ == refresh_coarse) {

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
