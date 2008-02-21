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

/**
 */

Array::Array(int  n0, int  n1, int  n2, int n3)
  : N_(0),
    a_(0)
{
  this->allocate_(n0,n1,n2,n3);
}

//----------------------------------------------------------------------

/**
 */

Array::~Array()
{
  this->deallocate_();
}

//----------------------------------------------------------------------

/**
 */
 
void Array::copy (const Array &array)
{
  this->deallocate_();
  this->allocate_(array.n_[0],array.n_[1],array.n_[2],array.n_[3]);
  this->copy_(array.values());
}

//----------------------------------------------------------------------

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

/**
 */

int Array::length () const
{
  return N_;
}


//----------------------------------------------------------------------

/**
 */

Scalar * Array::values () const
{
  return a_;
}

//----------------------------------------------------------------------

/**
 */

Scalar & Array::operator () (int i0, int i1, int i2, int i3)
{
  return a_[i0 + n_[0]*(i1 + n_[1]*(i2 +n_[2]*i3))];
}

//======================================================================


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
