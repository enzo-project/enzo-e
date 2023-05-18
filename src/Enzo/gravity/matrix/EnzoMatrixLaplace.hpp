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
  EnzoMatrixLaplace (int order = 4) throw()
    : mx_(0),
      my_(0),
      mz_(0),
      nx_(0),
      ny_(0),
      nz_(0),
      hx_(0.0),
      hy_(0.0),
      hz_(0.0),
      order_(order)
  {}

  /// Destructor
  virtual ~EnzoMatrixLaplace() throw()
  {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMatrixLaplace);

  /// CHARM++ migration constructor
  EnzoMatrixLaplace(CkMigrateMessage *m)
    : Matrix(m),
      mx_(0),
      my_(0),
      mz_(0),
      nx_(0),
      ny_(0),
      nz_(0),
      hx_(0.0),
      hy_(0.0),
      hz_(0.0),
      order_(0)
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
    p | order_;
  }

  /// Set cell widths.  Required for lower-level methods that don't have
  /// access to the Block
  void set_cell_width (double hx, double hy, double hz)
  {
    hx_ = hx;
    hy_ = hy;
    hz_ = hz;
  }
  
public: // virtual functions

  /// Apply the matrix to a vector Y <-- A*X
  virtual void matvec (int id_y, int id_x, Block * block, int g0=1) throw();

  /// Low-level matvec, useful for non-Block arrays (e.g. Block-local
  /// multigrid).  Must call set_cell_width and set_dimensions first
  /// manually!
  virtual void matvec (precision_type precision,
		       void * y, void * x, int g0=1) throw();

  /// Extract the diagonal into the given field
  virtual void diagonal (int id_x, Block * block, int g0=1) throw();

  /// Whether the matrix is singular or not
  virtual bool is_singular() const throw()
  { return true; }

  /// How many ghost zones required for matvec
  virtual int ghost_depth() const throw()
  { return (order_ == 2) ? 1 : ( (order_ == 4) ? 2 : 3); }

protected: // functions

  void matvec_ (enzo_float * Y, enzo_float * X, int g0) const throw();

  void diagonal_ (enzo_float * X, int g0) const throw();

protected: // attributes

  int mx_, my_, mz_;
  int nx_, ny_, nz_;
  double hx_, hy_, hz_;
  /// Order of the operator, 2 or 4
  int order_;

};

#endif /* COMPUTE_MATRIX_LAPLACE_HPP */
