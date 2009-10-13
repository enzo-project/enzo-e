//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      array.cpp
 * @brief     Declarations of Array member functions
 * @author    James Bordner
 * @date      Thu Feb 21 16:02:08 PST 2008
 *
 * $Id$
 * 
 *********************************************************************
 */
 
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "error_exception.hpp"
#include "array.hpp" 

//----------------------------------------------------------------------

/// Create a new uninitialized Array object

/**
 */

Array::Array() throw()
  : nx_(0),
    ny_(0),
    nz_(0),
    a_(0),
    is_allocated_(true)

{
}

//----------------------------------------------------------------------

/// Create an initialized array

/**
 */

Array::Array(Scalar * a,int nx,int ny,int nz) throw()
{
  nx_ = nx;
  ny_ = ny;
  nz_ = nz;
  is_allocated_ = false;
  a_ = a;
}

//----------------------------------------------------------------------

/// Deallocate the array

/**
 */

Array::~Array() throw()
{
  if (is_allocated_) {
    this->deallocate_();
    is_allocated_ = false;
  }
}

//----------------------------------------------------------------------

/// Copy an array into this one, deallocating any existing data

/**
 */
 
Array::Array (const Array &array) throw()
{
  this->reallocate_(array.nx_,array.ny_,array.nz_);
  this->copy_(array.values());
}

//----------------------------------------------------------------------

/// Assign an array to this one, deallocating any existing data

/**
 */
 
Array & Array::operator = (const Array &array) throw()
{
  this->reallocate_(array.nx_,array.ny_,array.nz_);
  this->copy_(array.values());
  return *this;
}

//----------------------------------------------------------------------

/// Resize the array, deallocating any existing data

/**
 */

void Array::resize (int nx, int ny, int nz) throw()
{
  if (nx != nx_ || 
      ny != ny_ || 
      nz != nz_) {
    this->reallocate_(nx,ny,nz);
  }
}


//----------------------------------------------------------------------

/// Return the size of the array

/**
 */

void Array::size (int * nx, int * ny, int * nz) const throw()
{
  if (nx) *nx = nx_;
  if (ny) *ny = ny_;
  if (nz) *nz = nz_;
}


//----------------------------------------------------------------------

/// Return the total length of the array

/**
 */

int Array::length () const throw()
{
  return nx_*ny_*nz_;
}


//----------------------------------------------------------------------

Scalar * Array::values () const throw()

/// Return a pointer to the array values

/**
 */

{
  return a_;
}

//----------------------------------------------------------------------

/// Return the given array element

/**
 */

Scalar & Array::operator () (int i0, int i1, int i2) throw()
{
  return a_[i0 + nx_*(i1 + ny_*i2)];
}

//----------------------------------------------------------------------

/// Set all values to 0, or to the given value if supplied

/**
 */

void Array::clear(Scalar value) throw()
{
  int i;
  for (i=0; i<nx_*ny_*nz_; i++) a_[i] = value;
}

//======================================================================

/**
 */

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

/**
 */

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

/**
 */

void Array::reallocate_(int nx, int ny, int nz) throw()
{
  deallocate_();
  allocate_(nx,ny,nz);
}

//----------------------------------------------------------------------

/**
 */

void Array::copy_(Scalar * a) throw()
{
  for (int i=0; i<nx_*ny_*nz_; i++) a_[i] = a[i];
}

//----------------------------------------------------------------------
