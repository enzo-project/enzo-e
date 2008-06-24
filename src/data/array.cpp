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

#include "string.h"
#include "scalar.hpp"
#include "io_error_.hpp"
#include "array.hpp" 

//----------------------------------------------------------------------

/// Create a new uninitialized Array object

/**
 */

Array::Array()
  : N_(0),
    a_(0)
{
  n_[0] = 0;
  n_[1] = 0;
  n_[2] = 0;
}

//----------------------------------------------------------------------

/// Create a new initialized Array object

/**
 */

Array::Array(int  n0, int  n1, int  n2)
  : N_(0),
    a_(0)
{
  this->allocate_(n0,n1,n2);
}

//----------------------------------------------------------------------

/// Deallocate the array

/**
 */

Array::~Array()
{
  this->deallocate_();
}

//----------------------------------------------------------------------

/// Copy an array into this one, deallocating any existing data

/**
 */
 
void Array::copy (const Array &array)
{
  this->deallocate_();
  this->allocate_(array.n_[0],array.n_[1],array.n_[2]);
  this->copy_(array.values());
}

//----------------------------------------------------------------------

/// Resize the array, deallocating any existing data

/**
 */

void Array::resize (int n0, int n1, int n2)
{
  if (n0 != n_[0] || 
      n1 != n_[1] || 
      n2 != n_[2]) {
    this->deallocate_();
    this->allocate_(n0,n1,n2);
  }
}


//----------------------------------------------------------------------

/// Return the size of the array

/**
 */

void Array::size (int * n0, int * n1, int * n2) const
{
  if (n0) *n0 = n_[0];
  if (n1) *n1 = n_[1];
  if (n2) *n2 = n_[2];
}


//----------------------------------------------------------------------

/// Return the total length of the array

/**
 */

int Array::length () const
{
  return N_;
}


//----------------------------------------------------------------------

Scalar * Array::values () const

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

Scalar & Array::operator () (int i0, int i1, int i2)
{
  return a_[i0 + n_[0]*(i1 + n_[1]*i2)];
}

//----------------------------------------------------------------------

/// Set all values to 0, or to the given value if supplied

/**
 */

void Array::clear(Scalar value)
{
  int i;
  for (i=0; i<N_; i++) a_[i] = value;
}

//======================================================================

/**
 */

void Array::allocate_(int n0, int n1, int n2)
{
  n_[0] = n0;
  n_[1] = n1;
  n_[2] = n2;

  N_ = n0*n1*n2;

  strcpy (warning_message,"Reallocating without deallocating");
  if (a_) WARNING_MESSAGE("Array::allocate_");

  a_ = new Scalar [N_];
}

//----------------------------------------------------------------------

/**
 */

void Array::deallocate_()
{
  delete [] a_;
  a_ = 0;
}

//----------------------------------------------------------------------

/**
 */

void Array::copy_(Scalar * a)
{
  for (int i=0; i<N_; i++) a_[i] = a[i];
}

//----------------------------------------------------------------------
