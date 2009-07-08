#ifndef ARRAY_ARRAY_HPP
#define ARRAY_ARRAY_HPP

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
 * @file      array_array_.hpp
 * @brief     Declaration of the Array abstract base class
 * @author    James Bordner
 * @date      Wed Jul  8 16:01:10 PDT 2009
 * @bug       none
 * @note      Adding Array base class, and renaming old Array to ArraySerial
 *
 * $Id$
 * 
 *********************************************************************
 */

class Array {

/** 
 *********************************************************************
 *
 * @class     Array
 * @brief     Define the interface for a 1D,2D, or 3D array
 * @ingroup   Array
 *
 * DEPENDENCIES
 *
 * Parallel: For controlling parallelism in distributed or threaded Arrays
 * Disk:     For I/O
 *
 *********************************************************************
 */

private:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

public:

  /// Create a new uninitialized Array object
  Array() throw() {};
  /// Deallocate the array
  ~Array() throw() {};
  /// Copy an array into this one, deallocating any existing data
  Array(const Array &) throw() {};
  /// Copy an array into this one, deallocating any existing data
  Array & operator = (const Array &) throw() {};
  /// Resize the array, deallocating any existing data
  virtual void resize (int n0, int n1=1, int n2=1) throw() = 0 ;
  /// Return the size of the array
  virtual void size (int * n0, int * n1=0, int * n2=0) const throw() = 0;
  /// Return the total length of the array
  virtual int  length() const throw() = 0;
  /// Return a pointer to the array values
  virtual Scalar * values () const throw() = 0;
  /// Return the given array element
  virtual Scalar & operator() (int  i0, int  i1=0, int  i2=0) throw() = 0;
  /// Set all values to 0, or to the given value if supplied
  virtual void clear(Scalar value = 0.0) throw() = 0;

};

#endif /* ARRAY_ARRAY_HPP */
