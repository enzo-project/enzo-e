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

  /// Compute residual R <-- B - A*X
  void residual (int ir, int ib, int ix, Block * block, int g0=1) throw();
  

public: // virtual functions

  /// Apply the matrix to a vector Y <-- A*X
  virtual void matvec (int iy, int ix, Block * block, int g0=1) throw() = 0;

  /// Extract the diagonal into the given field
  virtual void diagonal (int ix, Block * block, int g0=1) throw() = 0;

protected: // functions

  template<class T>
  void residual_ (T * ir, T * ib,
		 int mx, int my, int mz,
		 int ig0) throw();

};

#endif /* COMPUTE_MATRIX_HPP */
