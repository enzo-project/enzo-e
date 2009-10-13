#ifndef ARRAY_HPP
#define ARRAY_HPP

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

#include "cello.h"

/** 
 *********************************************************************
 *
 * @file      array.hpp
 * @brief     Declaration of the Array abstract base class
 * @author    James Bordner
 * @date      Wed Jul  8 16:01:10 PDT 2009
 * @bug       none
 * @note      Adding Array base class, and renaming old Array to ArraySerial
 *
 * $Id: array.hpp 715 2009-07-08 23:48:09Z bordner $
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

//   //--------------------------------------------------
//   // CONSTRUCTORS AND DESTRUCTORS
//   //--------------------------------------------------

//   /// Create a new uninitialized Array object
//   Array() throw() {};

//   /// Deallocate the array
//   ~Array() throw() {};

//   /// Copy an array into this one, deallocating any existing data
//   Array(const Array &) throw() {};

//   //--------------------------------------------------
//   // COPY, CLEAR
//   //--------------------------------------------------

//   /// Copy an array into this one, deallocating any existing data
//   Array & operator = (const Array &) throw() {};

//   /// Clear the Array to the given value
//   virtual void clear(Scalar value = 0.0) throw() = 0;
    	
//   //--------------------------------------------------
//   // ALLOCATION AND DEALLOCATION
//   //--------------------------------------------------

//   ///	Allocate storage for the Array
//   void allocate();

//   ///	Allocate storage for the Array
//   void deallocate();

//   /// 	Return whether storage is allocated for the Array
//   bool is_allocated();

//   /// Assert whether values should be deallocated when Array is deleted
//   void set_allocated();

//   //--------------------------------------------------
//   // ARRAY LAYOUT AND ELEMENT ACCESS
//   //--------------------------------------------------
    	
//   /// Set the size (and dimension) of the Array
//   void set_layout(Layout layout);

//   /// Get the size (and dimension) of the Array
//   Layout get_layout();

//   /// Set the Array elements to an existing array
//   void set_array();

//   ///	Return the Array elements as an array
//   Scalar get_array();

//   //--------------------------------------------------
//   // RESHAPE, SPLIT, MERGE
//   //--------------------------------------------------

//   /// 	Shrink the Array by some number of zones along each axis
//   void shrink();

//   /// 	Enlarge the Array by some number of zones along each axis
//   void grow();

//   /// Split an Array into two at some point along some axis
//   void split();

//   /// 	Merge two Arrays into one along some axis
//   void merge();

//   /// Resize the array, deallocating any existing data
//   virtual void resize (int n0, int n1=1, int n2=1) throw() = 0 ;

//   /// Return the size of the array
//   virtual void get_size (int * n0, int * n1=0, int * n2=0) const throw() = 0;

//   //--------------------------------------------------
//   // INPUT / OUTPUT
//   //--------------------------------------------------

//   ///	Write the Array to disk
//   write();

//   ///	Read the Array from disk
//   read();
    	
//   //--------------------------------------------------
//   // ITERATOR
//   //--------------------------------------------------

//   ///	Return an iterator over all Blocks of the Array
//   ItArrayBlock iterator();

//   // ELEMENT ACCESS

//   /// Return a pointer to the array values
//   virtual Scalar * values () const throw() = 0;
//   /// Return the given array element
//   virtual Scalar & operator() (int  i0, int  i1=0, int  i2=0) throw() = 0;

};   

#include "array_serial.hpp" /* Serial Array */
#include "block.hpp"

#endif /* ARRAY_HPP */
