// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     field_Array.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-02-20
/// @brief    Implementation of the Array class
 
#include "cello.hpp"

#include "error.hpp"
#include "field.hpp" 

//----------------------------------------------------------------------

Array::Array() throw()
  : nx_(0),
    ny_(0),
    nz_(0),
    a_(0),
    is_allocated_(true)

{
}

//----------------------------------------------------------------------

Array::Array(Scalar * a,int nx,int ny,int nz) throw()
{
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  is_allocated_ = false;
  a_ = a;
}

//----------------------------------------------------------------------

Array::~Array() throw()
{
  if (is_allocated_) {
    deallocate_();
    is_allocated_ = false;
  }
}

//----------------------------------------------------------------------

Array::Array (const Array &array) throw()
{
  reallocate_(array.nx_,array.ny_,array.nz_);
  copy_(array.values());
}

//----------------------------------------------------------------------

Array & Array::operator = (const Array &array) throw()
{
  reallocate_(array.nx_,array.ny_,array.nz_);
  copy_(array.values());
  return *this;
}

//----------------------------------------------------------------------

void Array::resize (int nx, int ny, int nz) throw()
{
  if (nx != nx_ || 
      ny != ny_ || 
      nz != nz_) {
    reallocate_(nx,ny,nz);
  }
}


//----------------------------------------------------------------------

void Array::size (int * nx, int * ny, int * nz) const throw()
{
  if (nx) *nx = nx_;
  if (ny) *ny = ny_;
  if (nz) *nz = nz_;
}


//----------------------------------------------------------------------

int Array::length () const throw()
{
  return nx_*ny_*nz_;
}


//----------------------------------------------------------------------

Scalar * Array::values () const throw()

{
  return a_;
}

//----------------------------------------------------------------------

Scalar & Array::operator () (int i0, int i1, int i2) throw()
{
  return a_[i0 + nx_*(i1 + ny_*i2)];
}

//----------------------------------------------------------------------

void Array::clear(Scalar value) throw()
{
  int i;
  for (i=0; i<nx_*ny_*nz_; i++) a_[i] = value;
}

//======================================================================

void Array::deallocate_() throw (ExceptionBadArrayDeallocation)
{
  if (is_allocated_) {
    delete [] a_;
  } else {
    throw ExceptionBadArrayDeallocation();
//     strcpy (warning_message,"Trying to deallocated non-allocated array");
//     WARNING_MESSAGE("Array::deallocate_");
  }
  a_ = 0;
}

//----------------------------------------------------------------------

void Array::allocate_(int nx, int ny, int nz) throw(ExceptionBadArrayAllocation)
{
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;

  if (a_ != NULL) {
    throw ExceptionBadArrayAllocation();
//     strcpy (warning_message,"Reallocating without deallocating");
//     WARNING_MESSAGE("Array::allocate_");
  }

  a_ = new Scalar [nx_*ny_*nz_];
}

//----------------------------------------------------------------------

void Array::reallocate_(int nx, int ny, int nz) throw()
{
  deallocate_();
  allocate_(nx,ny,nz);
}

//----------------------------------------------------------------------

void Array::copy_(Scalar * a) throw()
{
  for (int i=0; i<nx_*ny_*nz_; i++) a_[i] = a[i];
}

//----------------------------------------------------------------------
