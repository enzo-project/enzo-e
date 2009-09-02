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

/// Create a new uninitialized ArraySerial object

/**
 */

ArraySerial::ArraySerial() throw()
  : nx_(0),
    ny_(0),
    nz_(0),
    is_allocated_(true),
    a_(0)
{
}

//----------------------------------------------------------------------

/// Create an initialized array

/**
 */

ArraySerial::ArraySerial(Scalar * a,int nx,int ny,int nz) throw()
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

ArraySerial::~ArraySerial() throw()
{
  this->deallocate_();
}

//----------------------------------------------------------------------

/// Copy an array into this one, deallocating any existing data

/**
 */
 
ArraySerial::ArraySerial (const ArraySerial &array) throw()
{
  this->reallocate_(array.nx_,array.ny_,array.nz_);
  this->copy_(array.values());
}

//----------------------------------------------------------------------

/// Assign an array to this one, deallocating any existing data

/**
 */
 
ArraySerial & ArraySerial::operator = (const ArraySerial &array) throw()
{
  this->reallocate_(array.nx_,array.ny_,array.nz_);
  this->copy_(array.values());
  return *this;
}

//----------------------------------------------------------------------

/// Resize the array, deallocating any existing data

/**
 */

void ArraySerial::resize (int n0, int n1, int n2) throw()
{
  if (n0 != nx_ || 
      n1 != ny_ || 
      n2 != nz_) {
    this->reallocate_(n0,n1,n2);
  }
}


//----------------------------------------------------------------------

/// Return the size of the array

/**
 */

void ArraySerial::size (int * n0, int * n1, int * n2) const throw()
{
  if (n0) *n0 = nx_;
  if (n1) *n1 = ny_;
  if (n2) *n2 = nz_;
}


//----------------------------------------------------------------------

/// Return the total length of the array

/**
 */

int ArraySerial::length () const throw()
{
  return nx_*ny_*nz_;
}


//----------------------------------------------------------------------

Scalar * ArraySerial::values () const throw()

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

Scalar & ArraySerial::operator () (int i0, int i1, int i2) throw()
{
  return a_[i0 + nx_*(i1 + ny_*i2)];
}

//----------------------------------------------------------------------

/// Set all values to 0, or to the given value if supplied

/**
 */

void ArraySerial::clear(Scalar value) throw()
{
  int i;
  for (i=0; i<nx_*ny_*nz_; i++) a_[i] = value;
}

//======================================================================

/**
 */

void ArraySerial::deallocate_() throw (ExceptionBadArrayDeallocation)
{
  if (is_allocated_) {
    delete [] a_;
  } else {
    throw ExceptionBadArrayDeallocation();
//     strcpy (warning_message,"Trying to deallocated non-allocated array");
//     WARNING_MESSAGE("ArraySerial::deallocate_");
  }
  a_ = 0;
}

//----------------------------------------------------------------------

/**
 */

void ArraySerial::allocate_(int n0, int n1, int n2) throw(ExceptionBadArrayAllocation)
{
  nx_ = n0;
  ny_ = n1;
  nz_ = n2;

  if (a_ != NULL) {
    throw ExceptionBadArrayAllocation();
//     strcpy (warning_message,"Reallocating without deallocating");
//     WARNING_MESSAGE("ArraySerial::allocate_");
  }

  a_ = new Scalar [nx_*ny_*nz_];
}

//----------------------------------------------------------------------

/**
 */

void ArraySerial::reallocate_(int n0, int n1, int n2) throw()
{
  deallocate_();
  allocate_(n0,n1,n2);
}

//----------------------------------------------------------------------

/**
 */

void ArraySerial::copy_(Scalar * a) throw()
{
  for (int i=0; i<nx_*ny_*nz_; i++) a_[i] = a[i];
}

//----------------------------------------------------------------------
