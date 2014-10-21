// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityCg.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-21 17:25:40
/// @brief    [\ref Enzo] Declaration of EnzoMethodGravityCg
///
/// Congugate gradient (CG) method for solving for self-gravity on
/// fields.
///

#ifndef ENZO_ENZO_METHOD_GRAVITY_CG_HPP
#define ENZO_ENZO_METHOD_GRAVITY_CG_HPP

class EnzoMethodGravityCg : public Method {

  /// @class    EnzoMethodGravityCg
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve self-gravity
  /// using the CG method.  This is more applicable to smaller problems
  /// since the method doesn't scale as well as some other methods
  /// (FFT, MG, etc.) for larger problems.

public: // interface

  /// Create a new EnzoMethodGravityCg object
  EnzoMethodGravityCg(int iter_max, double res_tol);

  EnzoMethodGravityCg() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityCg);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityCg (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Solve for the gravitational potential
  virtual void compute( CommBlock * comm_block) throw();

protected: // methods

  template <class T>
  void compute_ (T * density,   int md3[3], int nd3[3],
		 T * potential, int mp3[3], int np3[3]) const throw();

protected: // attributes

  /// Maximum number of CG iterations
  int iter_max_;

  /// Convergence tolerance on the residual 
  double res_tol_;

};

#endif /* ENZO_ENZO_METHOD_GRAVITY_CG_HPP */
