// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityCg.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-21 17:25:40
/// @brief    [\ref Enzo] Declaration of EnzoMethodGravityCg
///
/// Bicongugate gradient stabalized method (Cg) for solving for
/// self-gravity on field data.

#ifndef ENZO_ENZO_METHOD_GRAVITY_CG_HPP
#define ENZO_ENZO_METHOD_GRAVITY_CG_HPP

class EnzoMethodGravityCg : public Method {

  /// @class    EnzoMethodGravityCg
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve self-gravity
  /// using the Cg method.  This is more applicable to smaller problems
  /// since the method doesn't scale as well as some other methods
  /// (FFT, MG, etc.) for larger problems.

public: // interface

  /// Create a new EnzoMethodGravityCg object
  EnzoMethodGravityCg(const FieldDescr * field_descr,
		      int iter_max, double res_tol);

  EnzoMethodGravityCg() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityCg);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityCg (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Solve for the gravitational potential
  virtual void compute( CommBlock * comm_block) throw();

  /// Continuation after global reduction
  void cg_loop_1(double rr) throw();

protected: // methods

  template <class T>
  void compute_ (EnzoBlock * enzo_block) throw();

  template <class T>
  void cg_begin_1_() throw();

  template <class T>
  void cg_exit_() throw();

protected: // attributes

  /// Maximum number of Cg iterations
  int iter_max_;

  /// Convergence tolerance on the residual 
  double res_tol_;

  /// Precision
  int precision_;

  /// Density and potential field id's

  int idensity_;
  int ipotential_;

  /// CG temporary field id's
  int ib_;
  int ix_;
  int ir_;
  int ip_;
  int iw_;
  int iap_;

  /// vector attributes
  int nx_,ny_,nz_;
  int mx_,my_,mz_;
  int gx_,gy_,gz_;

  /// CG variables
  int iter_;
  double alpha_;
  double pap_;
  double rr_;
  double rr_new_;

  
};

#endif /* ENZO_ENZO_METHOD_GRAVITY_CG_HPP */
