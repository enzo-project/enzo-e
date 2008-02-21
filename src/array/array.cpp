// $Id$
/**
 * @file    array.cpp
 * @brief   Member functions for the Array class
 * @author  James Bordner 
 * @version 1.0
 *
 * Attributes: see array.hpp
 * Operations: see array.hpp
 *
 */
// $Log$
 
#include <stdio.h>
#include <assert.h>

#include "scalar.hpp"
#include "error.hpp"
#include "array.hpp" 

//----------------------------------------------------------------------

/// Create an uninitialized Array
 
/**
 */

Array::Array()
  : N_(0),
    a_(0)
{
  n_[0] = 0;
  n_[1] = 0;
  n_[2] = 0;
  n_[3] = 0;
}

//----------------------------------------------------------------------

/// Create an initialized Array
 
/**
 */

Array::Array(int  n0, int  n1, int  n2, int n3)
  : N_(0),
    a_(0)
{
  this->allocate_(n0,n1,n2,n3);
}

//----------------------------------------------------------------------

/// Deallocate the Array
 
/**
 */
 
Array::~Array()
{
  this->deallocate_();
}

//----------------------------------------------------------------------

/// Copy an existing array

/**
 */
 
void Array::copy (const Array &array)
{
  this->deallocate_();
  this->allocate_(array.n_[0],array.n_[1],array.n_[2],array.n_[3]);
  this->copy_(array.values());
}

//----------------------------------------------------------------------

/// Set the dimension and extents of the array and allocate

/**
 */

void Array::resize (int n0, int n1, int n2, int n3)
{
  if (n0 != n_[0] || 
      n1 != n_[1] || 
      n2 != n_[2] || 
      n3 != n_[3]) {
    this->deallocate_();
    this->allocate_(n0,n1,n2,n3);
  }
}


//----------------------------------------------------------------------

/// Get the shape of the array

/**
 */

void Array::size (int * n0, int * n1, int * n2, int * n3) const
{
  if (n0) *n0 = n_[0];
  if (n1) *n1 = n_[1];
  if (n2) *n2 = n_[2];
  if (n3) *n3 = n_[3];
}


//----------------------------------------------------------------------

/// Return the length of the array, or 0 if uninitialized

/**
 */

int Array::length () const
{
  return N_;
}


//----------------------------------------------------------------------

/// Return a pointer to the array values (dangling pointer hazzard!)

/**
 */

Scalar * Array::values () const
{
  return a_;
}

//----------------------------------------------------------------------

/// Return the (i0,i1,i2,i3) element

/**
 */

Scalar & Array::operator () (int i0, int i1, int i2, int i3)
{
  return a_[i0 + n_[0]*(i1 + n_[1]*(i2 +n_[2]*i3))];
}

//======================================================================


/// Allocate array and initialize Array attributes

/**
 */

void Array::allocate_(int n0, int n1, int n2, int n3)
{
  n_[0] = n0;
  n_[1] = n1;
  n_[2] = n2;
  n_[3] = n3;

  N_ = n0*n1*n2*n3;

  if (a_) WARNING("Array::allocate_() reallocating without deallocating");

  a_ = new Scalar [N_];
}

//----------------------------------------------------------------------

void Array::deallocate_()
{
  delete [] a_;
  a_ = 0;
}

//----------------------------------------------------------------------

void Array::copy_(Scalar * a)
{
  for (int i=0; i<N_; i++) a_[i] = a[i];
}

//----------------------------------------------------------------------
