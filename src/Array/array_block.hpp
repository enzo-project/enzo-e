//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//

#ifndef ARRAY_BLOCK_HPP
#define ARRAY_BLOCK_HPP

/// Interface class between Array and low-level (C/fortran) routines.

/// Defines up to a 4-D fortran-like array for storing 1 or more 3D
/// arrays.  Axes can be permuted, including the index selecting the
/// array for storing interleaved arrays.
///
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Mon Oct 12 14:38:21 PDT 2009
/// @ingroup   Array
/// @note      

class Block {

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

public:

  /// Create a new uninitialized Block object
  Block() throw() {};

  /// Create a new initialized Block object
  Block(Scalar * values, 
	int * permute,
	int ndx,  int ndy,  int ndz,  int nda=1,
	int nx=0, int ny=0, int nz=0, int na=0)
    throw() :
    values_(values)
  { 
    if (permute) {
      p_[permute[0]] = 0; // default 1 = x
      p_[permute[1]] = 1; // default 2 = y
      p_[permute[2]] = 2; // default 3 = z
      p_[permute[3]] = 3; // default 4 = a
    } else {
      p_[0] = 0;
      p_[1] = 1;
      p_[2] = 2;
      p_[3] = 3;
    }

    nd_[p_[0]] = ndx;
    nd_[p_[1]] = ndy;
    nd_[p_[2]] = ndz;
    nd_[p_[3]] = nda;

    n_[p_[0]] = nx ? nx : ndx;
    n_[p_[1]] = ny ? ny : ndy;
    n_[p_[2]] = nz ? nz : ndz;
    n_[p_[3]] = na ? na : nda;

    m_[p_[0]] = n_[p_[3]];
    m_[p_[1]] = n_[p_[0]]*m_[p_[0]];
    m_[p_[2]] = n_[p_[1]]*m_[p_[1]];
    m_[p_[3]] = n_[p_[2]]*m_[p_[2]];

  }

  /// Return a pointer to the ith array
  Scalar * values(int i=0) const throw()
  {
    return & values_[i*m_[3]];
  }

  /// Return allocated block dimensions
  void get_dim 
  (int * ndx, 
   int * ndy = NULL,
   int * ndz = NULL) const throw()
  { 
    if (ndx) *ndx = nd_[p_[0]];
    if (ndy) *ndy = nd_[p_[1]];
    if (ndz) *ndz = nd_[p_[2]];
  }

  /// Return the array size
  void get_size 
  (int * nx,
   int * ny = NULL,
   int * nz = NULL, 
   int  *na = NULL) const throw()
  {
    if (nx) *nx = n_[p_[0]];
    if (ny) *ny = n_[p_[1]];
    if (nz) *nz = n_[p_[2]];
    if (na) *na = n_[p_[3]];
  }

  /// Return increments for loop index calculations
  void get_inc 
  (int * mx, 
   int * my = NULL,
   int * mz = NULL,
   int * ma = NULL) const throw () 	
  {
    if (mx) *mx = m_[p_[0]];
    if (my) *my = m_[p_[1]];
    if (mz) *mz = m_[p_[2]];
    if (ma) *ma = m_[p_[3]];
  }

private:

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

#endif /* ARRAY_BLOCK_HPP */
