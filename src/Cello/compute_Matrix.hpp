// See LICENSE_CELLO file for license and copyright information

/// @file     compute_Matrix.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Compute] Declaration for the Matrix class

#ifndef COMPUTE_MATRIX_HPP
#define COMPUTE_MATRIX_HPP

class Matrix : public PUP::able 
{
  /// @class    Matrix
  /// @ingroup  Compute
  /// @brief    [\ref Compute] Interface to an application compute / analysis / visualization function.

public: // interface

  /// Create a new Matrix
  Matrix () throw()
  {}

  /// Destructor
  virtual ~Matrix() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_abstract(Matrix);

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);
  }

public: // virtual functions

  /// Apply the matrix to a vector Y <-- A*X

  virtual void matvec (void * Y, const void * X) const throw() = 0;

protected: // functions

};

#endif /* COMPUTE_MATRIX_HPP */
