#ifndef BLOCK_HPP
#define BLOCK_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

#include "cello.h"

/** 
 *********************************************************************
 *
 * @file      block.hpp
 * @brief     Declaration of the Block class
 * @author    James Bordner
 * @date      Mon Oct 12 14:38:21 PDT 2009
 * @bug       none
 * @note      
 *
 * $Id$
 * 
 *********************************************************************
 */

class Block {

/** 
 *********************************************************************
 *
 * @class     Block
 * @brief     Interface class between Array and low-level (C/fortran) routines
 * @ingroup   Array
 *
 * DEPENDENCIES
 *
 *
 *********************************************************************
 */

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

public:

  /// Create a new uninitialized Block object
  Block() throw() {};
  /// Create a new uninitialized Block object
  Block(Scalar * array, int ndx,int nx,int ndy=1,int ny=0, int ndz=1,int nz=0)
    throw() :
    array_(array),
    ndx_(ndx),
    ndy_(ndy),
    ndz_(ndz),
    nx_(nx),
    ny_(ny),
    nz_(nz)
  { }
  /// Return the array dimensions
  void nd (int * ndx, int * ndy=0, int * ndz=0) const throw()
  { 
    if (ndx) *ndx = ndx_;
    if (ndy) *ndy = ndy_;
    if (ndz) *ndz = ndz_;
  }
  /// Return the array size
  void n (int * nx, int * ny=0, int * nz=0) const throw()
  {
    if (nx) *nx = nx_;
    if (ny) *ny = ny_;
    if (nz) *nz = nz_;
  }

  /// Return the total length of the block
  Scalar * array() const throw()
  {
    return array_;
  }

private:

  /// Pointer to the 0 element of the array
  Scalar * array_;

  /// Dimensions of the array
  int ndx_,ndy_,ndz_;

  /// Size of the array
  int nx_,ny_,nz_;

};   

#endif /* BLOCK_HPP */
