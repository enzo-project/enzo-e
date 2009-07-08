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

#include "error.hpp"
#include "array.hpp" 

//----------------------------------------------------------------------

/// Create a new uninitialized ArraySerial object

/**
 */

ArraySerial::ArraySerial() throw()
  : N_(0),
    a_(0)
{
  n_[0] = 0;
  n_[1] = 0;
  n_[2] = 0;
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
  this->deallocate_();
  this->allocate_(array.n_[0],array.n_[1],array.n_[2]);
  this->copy_(array.values());
}

//----------------------------------------------------------------------

/// Assign an array to this one, deallocating any existing data

/**
 */
 
ArraySerial & ArraySerial::operator = (const ArraySerial &array) throw()
{
  this->deallocate_();
  this->allocate_(array.n_[0],array.n_[1],array.n_[2]);
  this->copy_(array.values());
  return *this;
}

//----------------------------------------------------------------------

/// Resize the array, deallocating any existing data

/**
 */

void ArraySerial::resize (int n0, int n1, int n2) throw()
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

void ArraySerial::size (int * n0, int * n1, int * n2) const throw()
{
  if (n0) *n0 = n_[0];
  if (n1) *n1 = n_[1];
  if (n2) *n2 = n_[2];
}


//----------------------------------------------------------------------

/// Return the total length of the array

/**
 */

int ArraySerial::length () const throw()
{
  return N_;
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
  return a_[i0 + n_[0]*(i1 + n_[1]*i2)];
}

//----------------------------------------------------------------------

/// Set all values to 0, or to the given value if supplied

/**
 */

void ArraySerial::clear(Scalar value) throw()
{
  int i;
  for (i=0; i<N_; i++) a_[i] = value;
}

//======================================================================

/**
 */

void ArraySerial::allocate_(int n0, int n1, int n2) throw()
{
  n_[0] = n0;
  n_[1] = n1;
  n_[2] = n2;

  N_ = n0*n1*n2;

  strcpy (warning_message,"Reallocating without deallocating");
  if (a_) WARNING_MESSAGE("ArraySerial::allocate_");

  a_ = new Scalar [N_];
}

//----------------------------------------------------------------------

/**
 */

void ArraySerial::deallocate_() throw()
{
  delete [] a_;
  a_ = 0;
}

//----------------------------------------------------------------------

/**
 */

void ArraySerial::copy_(Scalar * a) throw()
{
  for (int i=0; i<N_; i++) a_[i] = a[i];
}

//----------------------------------------------------------------------
