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

class CProxy_ArrayMethodGravityCg;

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
		      int iter_max, double res_tol,
		      bool is_singular);

  EnzoMethodGravityCg() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityCg);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityCg (CkMigrateMessage *m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) ;
  
  /// Solve for the gravitational potential
  virtual void compute( Block * block) throw();

  /// Continuation after global reduction
  template <class T>
  void cg_loop_2(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  template <class T>
  void cg_loop_4(EnzoBlock * enzo_block) throw();

  /// Continuation after global reduction
  template <class T>
  void cg_loop_6(EnzoBlock * enzo_block) throw();

  template <class T>
  void cg_end (EnzoBlock * enzo_block, int retval) throw();

  /// Set rr_ by EnzoBlock after reduction
  void set_rr(long double rr) throw()  { rr_ = rr;  }

  /// Set bs_ by EnzoBlock after reduction
  void set_bs(long double bs) throw()  { bs_ = bs;  }

  /// Set bc_ by EnzoBlock after reduction
  void set_bc(long double bc) throw()  { bc_ = bc;  }

  /// Set pap_ by EnzoBlock after reduction
  void set_pap(double pap) throw()  { pap_ = pap; }

  /// Set iter_ by EnzoBlock after reduction
  void set_iter(int iter) throw()  { iter_ = iter; }

  /// Set rr_new_ by EnzoBlock after reduction
  void set_rr_new(double rr_new) throw()  { rr_new_ = rr_new; }

protected: // methods

  template <class T>
  void compute_ (EnzoBlock * enzo_block) throw();

  template <class T>
  void cg_begin_1_() throw();

  void cg_exit_() throw();

  /// Perform vector copy X <- Y
  template <class T>
  void copy_ (T * X, const T * Y) const throw();

  /// Compute local contribution to inner-product X*Y
  template <class T>
  long double dot_ (const T * X, const T * Y) const throw();

  template <class T>
  void zaxpy_ (T * Z, double a, const T * X, const T * Y) const throw();
  
  /// Compute local sum of vector elements X_i
  template <class T>
  long double sum_ (const T * X) const throw();

  /// scale the vector by the given scalar Y = a*X
  template <class T>
  void scale_ (T * Y, T a, const T * X) const throw();

  /// return the number of elements of the vector X
  template <class T>
  T count_ (T * X) const throw();
  
  /// Shift the vector X by a scalar multiple of Y
  /// NOTE includes ghost zones since performed after ghost refresh
  template <class T>
  void shift_ (T * X, const T a, const T * Y) const throw();

  /// Compute local matrix-vector product Y = A*X
  template <class T>
  void matvec_ (T * Y, const T * X) const throw();

  /// Set whether current Block is a leaf--if not don't touch data
  void set_leaf(Block * block) throw()
  { is_leaf_ = block->is_leaf(); }

protected: // attributes

  /// Whether you need to subtract of the nullspace of A from b, e.g. fully
  /// periodic or Neumann problems
  bool is_singular_;

  /// Whether current block is a leaf
  bool is_leaf_;

  /// Maximum number of Cg iterations
  int iter_max_;

  /// Convergence tolerance on the residual reduction rr_ / rr0_
  double res_tol_;

  /// Initial residual
  double rr0_;

  /// Mesh spacing for current Block
  double hx_, hy_, hz_;

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

  /// Current CG iteration
  int iter_;

  double alpha_;

  /// inner product P*(A*P)
  double pap_;

  /// inner product R*R
  long double rr_;

  /// newest inner product R*R
  double rr_new_;

  /// sum of elements B(i) for singular systems
  long double bs_;

  /// count of elements B(i) for singular systems
  long double bc_;

  /// block counter for operations that must be done once per processor, 
  /// e.g. increment iter_
  int block_count_;
  
};

//----------------------------------------------------------------------

class EnzoArrayMethodGravityCg : public CBase_EnzoArrayMethodGravityCg {

public:

  EnzoArrayMethodGravityCg () {};
  EnzoArrayMethodGravityCg (CkMigrateMessage*) {}; 

  /// EnzoMethodGravityCg synchronization entry method: DOT(R,R)
  template <class T>
  void r_cg_loop_1 (CkReductionMsg * msg) ;
  /// EnzoMethodGravityCg synchronization entry method: MATVEC (A,P)
  template <class T>
  void r_cg_loop_2 (CkReductionMsg * msg) ;
  /// EnzoMethodGravityCg synchronization entry method: DOT(P,AP)
  template <class T>
  void r_cg_loop_3 (CkReductionMsg * msg) ;
  /// EnzoMethodGravityCg synchronization entry method: DOT(R,R)
  template <class T>
  void r_cg_loop_4 (CkReductionMsg * msg) ;

};

#endif /* ENZO_ENZO_METHOD_GRAVITY_CG_HPP */
