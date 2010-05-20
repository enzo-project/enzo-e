// $Id: field_FieldBlock.cpp 1388 2010-04-20 23:57:46Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     field_FieldBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed May 19 18:17:50 PDT 2010
/// @brief    Implementation of the FieldBlock class

#include "field.hpp"

FieldBlock::FieldBlock() throw ()
  : values_(0)
{
  for (int i=0; i<4; i++) {
    nd_[i] = 0;
    n_[i] = 0;
    m_[i] = 0;
    p_[i] = 0;
  }
}

//----------------------------------------------------------------------

  /// Create a new initialized FieldBlock object
FieldBlock::FieldBlock(Scalar * values, 
		       int *    permute,
		       int      ndx,  
		       int      ndy,
		       int      ndz,
		       int      nda,
		       int      nx,
		       int      ny,
		       int      nz,
		       int      na) throw() 
  :    values_(values)
{ 
  if (permute) {
    p_[permute[0]] = 0; // default 1 = x
    p_[permute[1]] = 1; // default 2 = y
    p_[permute[2]] = 2; // default 3 = z
    p_[permute[3]] = 3; // default 4 = a
  } else {
    p_[0] = 0;
    p_[1] = 1;
    p_[2] = 2;
    p_[3] = 3;
  }

  nd_[p_[0]] = ndx;
  nd_[p_[1]] = ndy;
  nd_[p_[2]] = ndz;
  nd_[p_[3]] = nda;

  n_[p_[0]] = nx ? nx : ndx;
  n_[p_[1]] = ny ? ny : ndy;
  n_[p_[2]] = nz ? nz : ndz;
  n_[p_[3]] = na ? na : nda;

  m_[p_[0]] = n_[p_[3]];
  m_[p_[1]] = n_[p_[0]]*m_[p_[0]];
  m_[p_[2]] = n_[p_[1]]*m_[p_[1]];
  m_[p_[3]] = n_[p_[2]]*m_[p_[2]];

}

//----------------------------------------------------------------------

FieldBlock::~FieldBlock() throw ()
{
}

//----------------------------------------------------------------------

FieldBlock::FieldBlock(const FieldBlock & classname) throw ()
/// @param     classname  Object being copied
{
}

//----------------------------------------------------------------------

FieldBlock & FieldBlock::operator= (const FieldBlock & classname) throw ()
/// @param     classname  Source object of the assignment
/// @return    The target assigned object
{
  return *this;
}

//======================================================================

Scalar * FieldBlock::values(int i) const throw()
{
  return & values_[i*m_[3]];
}

//----------------------------------------------------------------------

void FieldBlock::get_dim 
(int * ndx, 
 int * ndy,
 int * ndz) const throw()
{ 
  if (ndx) *ndx = nd_[p_[0]];
  if (ndy) *ndy = nd_[p_[1]];
  if (ndz) *ndz = nd_[p_[2]];
}

//----------------------------------------------------------------------

void FieldBlock::get_size 
(int * nx,
 int * ny,
 int * nz, 
 int  *na) const throw()
{
  if (nx) *nx = n_[p_[0]];
  if (ny) *ny = n_[p_[1]];
  if (nz) *nz = n_[p_[2]];
  if (na) *na = n_[p_[3]];
}

//----------------------------------------------------------------------

void FieldBlock::get_inc 
(int * mx, 
 int * my,
 int * mz,
 int * ma) const throw () 	
{
  if (mx) *mx = m_[p_[0]];
  if (my) *my = m_[p_[1]];
  if (mz) *mz = m_[p_[2]];
  if (ma) *ma = m_[p_[3]];
}

//======================================================================
