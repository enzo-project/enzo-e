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

enum return_enum {
  return_unknown,
  return_converged,
  return_error_max_iter_reached
};

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
  template <class T>
  void cg_loop_1(EnzoBlock * enzo_block,double rr) throw();

  /// Continuation after global reduction
  template <class T>
  void cg_loop_2(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  template <class T>
  void cg_loop_3(EnzoBlock * enzo_block,double pap) throw();

  /// Continuation after global reduction
  template <class T>
  void cg_loop_4(EnzoBlock * enzo_block,double rr) throw();

  void cg_end (int retval) throw();

  /// Return the precision of the fields, used for calling templated methods
  int precision() const throw()
  { return precision_; }

protected: // methods

  template <class T>
  void compute_ (EnzoBlock * enzo_block) throw();

  template <class T>
  void cg_begin_1_() throw();

  void cg_exit_() throw();

  /// Compute local contribution to inner-product X*Y
  template <class T>
  T dot_ (const T * X, const T * Y) const throw();

  template <class T>
  void zaxpy_ (T * Z, double a, const T * X, const T * Y) const throw();
  
  /// Compute local matrix-vector product Y = A*X
  template <class T>
  void matvec_ (T * Y, const T * X) const throw();
   

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
  int iap_;

  /// vector attributes
  int nx_,ny_,nz_;
  int mx_,my_,mz_;
  int gx_,gy_,gz_;

  /// CG variables
  int iter_;
  double alpha_;
  double rr_;
  double rr_new_;

  
};

#endif /* ENZO_ENZO_METHOD_GRAVITY_CG_HPP */
