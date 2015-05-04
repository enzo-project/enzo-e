// See LICENSE_CELLO file for license and copyright information

/// @file     compute_EnzoMatrixDiagonal.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2015-04-02
/// @brief    [\ref Compute] Declaration of the EnzoMatrixDiagonal class

#ifndef COMPUTE_MATRIX_DIAGONAL_HPP
#define COMPUTE_MATRIX_DIAGONAL_HPP

class EnzoMatrixDiagonal : public Matrix 
{
  /// @class    EnzoMatrixDiagonal
  /// @ingroup  Compute
  /// @brief    [\ref Compute] Interface to an application compute / analysis / visualization function.

public: // interface

  /// Create a new EnzoMatrixDiagonal
  EnzoMatrixDiagonal () throw()
  {}

  /// Destructor
  virtual ~EnzoMatrixDiagonal() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMatrixDiagonal);

  /// CHARM++ migration constructor
  EnzoMatrixDiagonal(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);
    p | m_;
    p | hx_;
    p | hy_;
    p | hz_;
  }

public: // virtual functions

  /// Apply the matrix to a vector Y <-- A*X
  virtual void matvec (int id_y, int id_x, Block * block) throw();

  /// Extract the diagonal into the given field
  virtual void diagonal (int id_x, Block * block) throw();

protected: // functions

  template <class T>
  void matvec_ (T * Y, T * X) const throw();

  template <class T>
  void diagonal_ (T * X) const throw();

protected: // attributes

  double hx_, hy_, hz_;
  int m_;
};

#endif /* COMPUTE_MATRIX_DIAGONAL_HPP */
