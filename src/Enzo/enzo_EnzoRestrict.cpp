// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRestrict.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-10
/// @brief    Implentation of Enzo's restriction operators

#include "enzo.hpp"

#include "cello.hpp"

EnzoRestrict::EnzoRestrict(std::string restrict_type) throw()
  : Restrict ()
{
}
//----------------------------------------------------------------------

void EnzoRestrict::pup (PUP::er &p)
{
  TRACEPUP;

  Restrict::pup(p);
}

//----------------------------------------------------------------------

int EnzoRestrict::apply 
( precision_type precision,
  void *       values_c, int nd3_c[3], int im3_c[3],  int n3_c[3],
  const void * values_f, int nd3_f[3], int im3_f[3],  int n3_f[3])
{
  switch (precision)  {

  case precision_single:

    return apply_( (float *)       values_c, nd3_c, im3_c, n3_c,
		   (const float *) values_f, nd3_f, im3_f, n3_f);

    break;

  case precision_double:

    return apply_( (double *)       values_c, nd3_c, im3_c, n3_c,
		   (const double *) values_f, nd3_f, im3_f, n3_f);

    break;

  default:

    ERROR1 ("EnzoRestrict::apply()",
	    "Unknown precision %d",
	    precision);
    return 0;
  }
}

//----------------------------------------------------------------------

template <class T>
int EnzoRestrict::apply_
( T *       values_c, int nd3_c[3], int im3_c[3],  int n3_c[3],
  const T * values_f, int nd3_f[3], int im3_f[3],  int n3_f[3])
{
  return 0;
  // int nd3_f[3];
  // int nd3_c[3];
  // field_block_c->size(&nd3_c[0],&nd3_c[1],&nd3_c[2]);
  // field_block_f->size(&nd3_f[0],&nd3_f[1],&nd3_f[2]);

  // int rank = (nd3_c[2] > 1) ? 3 : (nd3_c[1] >1 ) ? 2 : 1;

  // for (int index=0; index<field_descr->field_count(); index++) {

  //   enzo_float * values_f = (enzo_float *) field_block_f->values(index);
  //   enzo_float * values_c = (enzo_float *) field_block_c->values(index);

  //   int i3m_f[3];
  //   int i3m_c[3];
  //   field_descr->ghosts(index,&i3m_c[0],&i3m_c[1],&i3m_c[2]);
  //   field_descr->ghosts(index,&i3m_f[0],&i3m_f[1],&i3m_f[2]);

  //   int r3[3];
  //   int i3p_f[3];
  //   int i3p_c[3];
  //   for (int axis=0; axis<3; axis++) {
  //     r3[axis] = (rank > axis) ? 2 : 1;
  //     i3p_f[axis]=nd3_f[axis]-i3m_f[axis];
  //     i3p_c[axis]=nd3_c[axis]-i3m_c[axis];
  //   }

  //   enzo_float * work = 0;
  //   int positivity_flag = 2;
  //   int error;

  //   // shift fine values to those over child

  //   int icx = icx*(i3p_f[0] - i3m_f[0])/2;
  //   int icy = icy*(i3p_f[1] - i3m_f[1])/2;
  //   int icz = icz*(i3p_f[2] - i3m_f[2])/2;

  //   int ic = icx + nd3_f[0]*(icy + nd3_f[1]*icz);

  //   values_f = &values_f[ic];


  //   ASSERT1("EnzoRestrict::apply",
  // 	    "apply() returned error %d",
  // 	    error, ! error);
  // }
}  
