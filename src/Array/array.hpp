// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ARRAY_HPP
#define ARRAY_HPP

/// @file     array.hpp
/// @brief    Include file for the Array component
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul  8 16:01:10 PDT 2009
/// @todo     Split file into array.hpp module include file and
///           array_array.hpp class definition file

#include "cello.h"
#include "error_exception.hpp"

class Array {

  /// @class    Array
  /// @brief    Defines the interface for arrays
  /// @ingroup  Array
  /// @todo     Add commented-out functions at bottom

public: // interface

  /// Initialize an empty Array object
  Array() throw();

  /// Initiazize a new Array object
  Array (Scalar * values,
	 int      nx,
	 int      ny=1,
	 int      nz=1) throw();

  /// Deallocate the Array object
  ~Array() throw();

  /// Initialize the Array as a copy of another
  Array(const Array &) throw();

  /// Copy an array into this one, deallocating any existing data
  Array & operator = (const Array &) throw();

  /// Resize the array, deallocating any existing data
  virtual void resize (int n0, 
		       int n1=1,
		       int n2=1) throw();

  /// Return the size of the array
  virtual void size (int * n0, 
		     int * n1=0, 
		     int * n2=0) const throw();

  /// Return the total length of the array
  virtual int length() const throw();

  /// Return a pointer to the array values
  virtual Scalar * values () const throw();

  /// Return the given array element
  virtual Scalar & operator() (int  i0, 
			       int  i1=0, 
			       int  i2=0) throw();

  /// Set all values to 0, or to the given value if supplied
  virtual void clear(Scalar value = 0.0) throw();

private: // functions

  /// Deallocate array values
  void deallocate_() throw(ExceptionBadArrayDeallocation);

  /// Allocate array values
  void allocate_ (int n0, 
		  int n1=1, 
		  int n2=1) throw(ExceptionBadArrayAllocation);

  /// Re-allocate array values
  void reallocate_ (int n0, 
		    int n1=1, 
		    int n2=1) throw();

  /// Copy array values a[0:N-1] to ArraySerial
  void copy_ (Scalar * a) throw();

private: // attributes

  /// Shape of the array, right-padded with 1's
  int nx_,ny_,nz_;

  /// Array values stored in column-major ordering
  Scalar *a_;

  /// Whether the array values are allocated internally or externally
  bool is_allocated_;

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

#include "array_block.hpp"

#endif /* ARRAY_HPP */
