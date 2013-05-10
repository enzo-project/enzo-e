// See LICENSE_CELLO file for license and copyright information

/// @file     field_EnzoProlong.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of Enzo's prolongation

#include "enzo.hpp"

EnzoProlong::EnzoProlong(std::string prolong_type) throw()
  : Prolong ()
{
  if      (prolong_type == "ThirdOrderA")  interpolation_method_ = 0;
  else if (prolong_type == "SecondOrderA") interpolation_method_ = 1;
  else if (prolong_type == "SecondOrderB") interpolation_method_ = 2;
  else if (prolong_type == "SecondOrderC") interpolation_method_ = 3;
  else if (prolong_type == "FirstOrderA")  interpolation_method_ = 4;
  else {
    ERROR1("EnzoProlong::EnzoProlong",
	  "Unrecognized interpolation method %s",
	   prolong_type.c_str());
  }
}
//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void EnzoProlong::pup (PUP::er &p)
{
  TRACEPUP;

  Prolong::pup(p);

  p | interpolation_method_;
}

#endif /* CONFIG_USE_CHARM */

//----------------------------------------------------------------------

void EnzoProlong::apply 
(
 CommBlock        * comm_block_h, 
 const  CommBlock * comm_block_H, 
 const FieldDescr * field_descr,
 int icx, int icy, int icz)
{
  Block * block_h = comm_block_h->block();
  const Block * block_H = comm_block_H->block();
  FieldBlock * field_block_h = block_h->field_block();
  const FieldBlock * field_block_H = block_H->field_block();

  for (int index=0; index<field_descr->field_count(); index++) {

    int gx,gy,gz;
    field_descr->ghosts(index,&gx,&gy,&gz);

    int nx,ny,nz;
    field_block_h->size(&nx,&ny,&nz);

     //  subroutine interpolate(ndim, parent, pdims, pstart, pend, refine,
     // &                       grid, gdims, gstart, work, imethod,
     // &                       iposflag, ierror)
      

  }
  
}

//======================================================================

