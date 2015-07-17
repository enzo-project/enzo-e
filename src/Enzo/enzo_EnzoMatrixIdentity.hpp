// See LICENSE_CELLO file for license and copyright information

/// @file     compute_EnzoMatrixIdentity.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2015-04-02
/// @brief    [\ref Compute] Declaration of the EnzoMatrixIdentity class

#ifndef COMPUTE_MATRIX_IDENTITY_HPP
#define COMPUTE_MATRIX_IDENTITY_HPP

class EnzoMatrixIdentity : public Matrix 
{
  /// @class    EnzoMatrixIdentity
  /// @ingroup  Compute
  /// @brief    [\ref Compute] Interface to an application compute / analysis / visualization function.

public: // interface

  /// Create a new EnzoMatrixIdentity
  EnzoMatrixIdentity () throw()
  {}

  /// Destructor
  virtual ~EnzoMatrixIdentity() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMatrixIdentity);

  /// CHARM++ migration constructor
  EnzoMatrixIdentity(CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  { TRACEPUP;
    PUP::able::pup(p);
    p | mx_;
    p | my_;
    p | mz_;
  }

public: // virtual functions

  /// Apply the matrix to a vector Y <-- A*X
  virtual void matvec (int id_y, int id_x, Block * block, int g0 = 1) throw();

  /// Extract the diagonal into the given field
  virtual void diagonal (int id_x, Block * block, int g0 = 1) throw();

protected: // functions

  template <class T>
  void matvec_ (T * Y, T * X, int g0) const throw();

  template <class T>
  void diagonal_ (T * X, int g0) const throw();

protected: // attributes

  int mx_,my_,mz_;

};

#endif /* COMPUTE_MATRIX_IDENTITY_HPP */
