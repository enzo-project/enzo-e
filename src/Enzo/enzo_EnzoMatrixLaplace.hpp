// See LICENSE_CELLO file for license and copyright information

/// @file     compute_EnzoMatrixLaplace.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2015-04-02
/// @brief    [\ref Compute] Declaration of the EnzoMatrixLaplace class

#ifndef COMPUTE_MATRIX_LAPLACE_HPP
#define COMPUTE_MATRIX_LAPLACE_HPP

class EnzoMatrixLaplace : public Matrix 
{
  /// @class    EnzoMatrixLaplace
  /// @ingroup  Compute
  /// @brief    [\ref Compute] Interface to an application compute / analysis / visualization function.

public: // interface

  /// Create a new EnzoMatrixLaplace
  EnzoMatrixLaplace () throw()
  {}

  /// Destructor
  virtual ~EnzoMatrixLaplace() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMatrixLaplace);

  /// CHARM++ migration constructor
  EnzoMatrixLaplace(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);
    p | mx_;
    p | my_;
    p | mz_;
    p | nx_;
    p | ny_;
    p | nz_;
    p | gx_;
    p | gy_;
    p | gz_;
    p | hx_;
    p | hy_;
    p | hz_;
    p | rank_;
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

  int mx_, my_, mz_;
  int nx_, ny_, nz_;
  int gx_, gy_, gz_;
  double hx_, hy_, hz_;
  int rank_;

};

#endif /* COMPUTE_MATRIX_LAPLACE_HPP */
