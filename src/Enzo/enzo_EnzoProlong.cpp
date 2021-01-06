// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoProlong.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-05-09
/// @brief    Implentation of Enzo's prolongation

#include "enzo.hpp"

#include "cello.hpp"

EnzoProlong::EnzoProlong(std::string type,int positive) throw()
  : Prolong (),
    method_(-1),
    positive_(positive)
{
  if      (type == "3A")  method_ = 0;
  else if (type == "2A") method_ = 1;
  else if (type == "2B") method_ = 2;
  else if (type == "2C") method_ = 3;
  else if (type == "1A")  method_ = 4;
  else {
    ERROR1("EnzoProlong::EnzoProlong",
	  "Unrecognized interpolation method %s",
	   type.c_str());
  }
}
//----------------------------------------------------------------------

void EnzoProlong::pup (PUP::er &p)
{
  TRACEPUP;

  Prolong::pup(p);

  p | method_;
}

//----------------------------------------------------------------------

int EnzoProlong::apply 
( precision_type precision,
  void *       values_f, int m3_f[3], int o3_f[3], int n3_f[3],
  const void * values_c, int m3_c[3], int o3_c[3], int n3_c[3],
  bool accumulate)
{
  return apply_((enzo_float *)     values_f,m3_f,o3_f,n3_f,
                (const enzo_float*)values_c,m3_c,o3_c,n3_c,accumulate);
}

//----------------------------------------------------------------------

int EnzoProlong::apply_
( 
 enzo_float * values_f, int m3_f[3], int o3_f[3], int n3_f[3],
 const enzo_float * values_c, int m3_c[3], int o3_c[3], int n3_c[3],
 bool accumulate)
{

  INCOMPLETE("EnzoProlong::apply()");

  int rank = cello::rank();

  int r3[3] = {2,2,2};
  int error = 0;
  int size=(n3_f[0]-o3_f[0]+1)/2 + 1;
  if (rank >= 2) size*=(n3_f[1]-o3_f[1]+1)/2 + 1;
  if (rank >= 3) size*=(n3_f[2]-o3_f[2]+1)/2 + 1;

  enzo_float * work = new enzo_float[size];
  
  FORTRAN_NAME(interpolate)
    (&rank,
     (enzo_float*)values_c, m3_c, o3_c, n3_c, r3,
     values_f, m3_f, o3_f, work, &method_,
     &positive_, &error);

  //--------------------------------------------------

  ASSERT1("EnzoProlong::apply",
          "apply() returned error %d",
          error, (! error) );
  return 0;
}  
