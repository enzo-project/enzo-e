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
  return apply_( (enzo_float *)       values_c, nd3_c, im3_c, n3_c,
		 (const enzo_float *) values_f, nd3_f, im3_f, n3_f);
}

//----------------------------------------------------------------------

int EnzoRestrict::apply_
( enzo_float *       values_c, int nd3_c[3], int im3_c[3],  int n3_c[3],
  const enzo_float * values_f, int nd3_f[3], int im3_f[3],  int n3_f[3])
{
  return 0;
}  
