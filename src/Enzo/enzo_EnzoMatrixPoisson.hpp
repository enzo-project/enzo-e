// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMatrixPoisson.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-27 22:37:41
/// @brief    [\ref Enzo] Implementation of Enzo's MatrixPoisson functions

#ifndef ENZO_ENZO_MATRIX_POISSON_HPP
#define ENZO_ENZO_MATRIX_POISSON_HPP

class EnzoMatrixPoisson : public Matrix {

  /// @class    EnzoMatrixPoisson
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's MatrixPoisson functions

public: // interface

  /// Create a new EnzoMatrixPoisson object
  EnzoMatrixPoisson () {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMatrixPoisson);
  
  /// Charm++ PUP::able migration constructor
  EnzoMatrixPoisson (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) {};
  
  /// Perform the computation on the block
  virtual void matvec( void * Y, const void * X) const throw()
  { printf ("%s:%d matvec()\n",__FILE__,__LINE__); }

private: // functions

  template <typename T>
  void matvec_(T * Y, const T * X) throw()
  { };

private: // attributes

};

#endif /* ENZO_ENZO_COMPUTE_POISSON_HPP */
