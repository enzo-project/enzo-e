// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMatrixLaplace.hpp 
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
    : mx_(0),
      my_(0),
      mz_(0),
      nx_(0),
      ny_(0),
      nz_(0),
      hx_(0.0),
      hy_(0.0),
      hz_(0.0),
      rank_(0)
  {}

  /// Destructor
  virtual ~EnzoMatrixLaplace() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMatrixLaplace);

  /// CHARM++ migration constructor
  EnzoMatrixLaplace(CkMigrateMessage *m)
    : mx_(0),
      my_(0),
      mz_(0),
      nx_(0),
      ny_(0),
      nz_(0),
      hx_(0.0),
      hy_(0.0),
      hz_(0.0),
      rank_(0)
  { }

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
    p | hx_;
    p | hy_;
    p | hz_;
    p | rank_;
  }

public: // virtual functions

  /// Apply the matrix to a vector Y <-- A*X
  virtual void matvec (int id_y, int id_x, Block * block, int g0=1) throw();

  /// Extract the diagonal into the given field
  virtual void diagonal (int id_x, Block * block, int g0=1) throw();

  /// Whether the matrix is singular or not
  virtual bool is_singular() const throw()
  { return true; }

protected: // functions

  template <class T>
  void matvec_ (T * Y, T * X, int g0) const throw();

  template <class T>
  void diagonal_ (T * X, int g0) const throw();

protected: // attributes

  int mx_, my_, mz_;
  int nx_, ny_, nz_;
  double hx_, hy_, hz_;
  int rank_;

};

#endif /* COMPUTE_MATRIX_LAPLACE_HPP */
