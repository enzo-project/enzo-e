// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef FIELD_BLOCK_HPP
#define FIELD_BLOCK_HPP

/// @file     field_FieldBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Oct 12 14:38:21 PDT 2009
/// @todo     Add allocation and deallocation
/// @todo     Implement and test merge(),split()
/// @todo     Implement and test grow(),shrink()
/// @todo     Implement and test read(),write()
/// @brief    Fortran-style array class.

class FieldBlock {

  /// @class    FieldBlock
  /// @ingroup  Field
  /// @brief    Interface between field arrays and low-level (C/fortran) routines.
  ///
  /// Defines up to a 4-D fortran-like array for storing 1 or more 3D
  /// arrays.  Axes can be permuted, including the index selecting the
  /// array for storing interleaved arrays.

public: // interface

  /// Create a new uninitialized FieldBlock object
  FieldBlock() throw();

  /// Create a new initialized FieldBlock object
  FieldBlock(Scalar * values, 
	     int *    permute,
	     int      ndx,  
	     int      ndy,
	     int      ndz,
	     int      nda = 1,
	     int      nx  = 0,
	     int      ny  = 0,
	     int      nz  = 0,
	     int      na  = 0) throw();

  /// Deconstructor
  ~FieldBlock() throw();

  /// Copy constructor
  FieldBlock(const FieldBlock & classname) throw ();

  /// Assignment operator
  FieldBlock & operator= (const FieldBlock & classname) throw ();

  /// Return a pointer to the ith array
  Scalar * values(int i=0) const throw();

  /// Return allocated block dimensions
  void get_dim 
  (int * ndx, 
   int * ndy = 0,
   int * ndz = 0) const throw();

  /// Return the array size
  void get_size 
  (int * nx,
   int * ny = 0,
   int * nz = 0, 
   int  *na = 0) const throw();

  /// Return increments for loop index calculations
  void get_inc 
  (int * mx, 
   int * my = 0,
   int * mz = 0,
   int * ma = 0) const throw ();

  // /// Resize the array, deallocating any existing data
  // void resize (int n0, 
  // 		       int n1=1,
  // 		       int n2=1) throw();

  /// Return the total length of the array
  int length() const throw();

  /// Set all values to 0, or to the given value if supplied
  void clear(Scalar value = 0.0) throw();

  //  /// Shrink the Array by some number of zones along each axis.  Used
  //  /// for deallocating ghost or boundary zones
  //   void shrink();

  //  /// Enlarge the Array by some number of zones along each axis. Used
  //  /// for allocating ghost or boundary zones.
  //   void grow();

  //   /// Split Block into multiple blocks along some subset of axes.  Used for AMR.
  //   void split();

  //   ///  Merge Blocks into a single one along some subset of axes.  Used for AMR.
  //   void merge();

  //   /// Return a subarray of the given block.  Used for accessing ghost or boundary zones.
  //   void subarray();

  //   /// Write the Array to disk
  //   write();

  //   /// Read the Array from disk
  //   read();

private: // attributes

  /// Pointer to the first element of the array
  Scalar * values_;

  /// Dimensions of the array
  int nd_[4];

  /// Size of the array
  int n_[4];

  /// Increments along each dimension
  int m_[4];

  /// Permutation such that p_[0:3] = [ix, iy, iz, ia]
  int p_[4];

};   

#endif /* FIELD_BLOCK_HPP */
