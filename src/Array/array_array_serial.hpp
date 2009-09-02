#ifndef ARRAY_ARRAY_SERIAL_HPP
#define ARRAY_ARRAY_SERIAL_HPP

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
 * @file      array_array_serial_.hpp
 * @brief     Declaration of the Array class
 * @author    James Bordner
 * @date      Thu Feb 21 13:54:19 PST 2008
 * @bug       none
 *
 * $Id: array_array_.hpp 712 2009-07-08 22:55:25Z bordner $
 * 
 *********************************************************************
 */

class ArraySerial : public Array {

/** 
 *********************************************************************
 *
 * @class     ArraySerial
 * @brief     Encapsulate a serial fortran-style 1D,2D, or 3D array
 * @ingroup   Array
 *
 *********************************************************************
 */

public:

  /// Create a new uninitialized ArraySerial object
  ArraySerial() throw();
  /// Create an initialized array
  ArraySerial(Scalar * values,int nx,int ny=1,int nz=1) throw();
  /// Deallocate the array
  ~ArraySerial() throw();
  /// Copy an array into this one, deallocating any existing data
  ArraySerial(const ArraySerial &) throw();
  /// Copy an array into this one, deallocating any existing data
  ArraySerial & operator = (const ArraySerial &) throw();
  /// Resize the array, deallocating any existing data
  virtual void resize (int n0, int n1=1, int n2=1) throw();
  /// Return the size of the array
  virtual void size (int * n0, int * n1=0, int * n2=0) const throw();
  /// Return the total length of the array
  virtual int  length() const throw();
  /// Return a pointer to the array values
  virtual Scalar * values () const throw();
  /// Return the given array element
  virtual Scalar & operator() (int  i0, int  i1=0, int  i2=0) throw();
  /// Set all values to 0, or to the given value if supplied
  virtual void clear(Scalar value = 0.0) throw();

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

private:

  /// Deallocate array values
  void deallocate_() throw(ExceptionBadArrayDeallocation);

  /// Allocate array values
  void allocate_ (int n0, int n1=1, int n2=1) throw(ExceptionBadArrayAllocation);

  /// Allocate array values
  void reallocate_ (int n0, int n1=1, int n2=1) throw();

  /// Copy array values a[0:N-1] to ArraySerial
  void copy_ (Scalar * a) throw();

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

private:

  /// Shape of array, right-padded with 1's
  int     nx_,ny_,nz_;

  /// Array values stored in column-major ordering
  Scalar *a_;

  /// Whether the array values are allocated internally or externally
  bool is_allocated_;

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

};

#endif /* ARRAY_ARRAY_SERIAL_HPP */
