// See LICENSE_CELLO file for license and copyright information

/// @file     field_EnzoRestrict.cpp
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

#ifdef CONFIG_USE_CHARM

void EnzoRestrict::pup (PUP::er &p)
{
  TRACEPUP;

  Restrict::pup(p);

}

#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

void EnzoRestrict::apply 
(
 FieldBlock        * field_block_c, 
 const  FieldBlock * field_block_f, 
 const FieldDescr * field_descr,
 int icx, int icy, int icz)
{

  int nd3_f[3];
  int nd3_c[3];
  field_block_c->size(&nd3_c[0],&nd3_c[1],&nd3_c[2]);
  field_block_f->size(&nd3_f[0],&nd3_f[1],&nd3_f[2]);

  int rank = (nd3_c[2] > 1) ? 3 : (nd3_c[1] >1 ) ? 2 : 1;

  for (int index=0; index<field_descr->field_count(); index++) {

    enzo_float * values_f = (enzo_float *) field_block_f->field_values(index);
    enzo_float * values_c = (enzo_float *) field_block_c->field_values(index);

    int i3m_f[3];
    int i3m_c[3];
    field_descr->ghosts(index,&i3m_c[0],&i3m_c[1],&i3m_c[2]);
    field_descr->ghosts(index,&i3m_f[0],&i3m_f[1],&i3m_f[2]);

    int r3[3];
    int i3p_f[3];
    int i3p_c[3];
    for (int axis=0; axis<3; axis++) {
      r3[axis] = (rank > axis) ? 2 : 1;
      i3p_f[axis]=nd3_f[axis]-i3m_f[axis];
      i3p_c[axis]=nd3_c[axis]-i3m_c[axis];
    }

    enzo_float * work = 0;
    int positivity_flag = 2;
    int error;

    // shift fine values to those over child

    int icx = icx*(i3p_f[0] - i3m_f[0])/2;
    int icy = icy*(i3p_f[1] - i3m_f[1])/2;
    int icz = icz*(i3p_f[2] - i3m_f[2])/2;

    int ic = icx + nd3_f[0]*(icy + nd3_f[1]*icz);

    values_f = &values_f[ic];


    ASSERT1("EnzoRestrict::apply",
	    "apply() returned error %d",
	    error, ! error);
  }
}  
