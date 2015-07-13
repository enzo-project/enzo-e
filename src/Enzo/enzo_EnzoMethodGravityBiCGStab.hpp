// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityBiCGStab.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @date     2014-10-21 17:25:40
/// @brief    [\ref Enzo] Declaration of EnzoMethodGravityBiCGStab
///
/// Bicongugate gradient stabilized method (BiCGStab) for solving 
/// for self-gravity on field data.

#ifndef ENZO_ENZO_METHOD_GRAVITY_BICGSTAB_HPP
#define ENZO_ENZO_METHOD_GRAVITY_BICGSTAB_HPP

class EnzoMethodGravityBiCGStab : public Method {

  /// @class    EnzoMethodGravityBiCGStab
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve self-gravity
  /// using the BiCGStab method.  Alone, this is more applicable to 
  /// smaller problems since the method doesn't scale as well as 
  /// some other methods (FFT, MG, etc.) for larger problems.  
  /// Alternately, a more scalable method may be combined as a 
  /// preconditioner for a robust and scalable overall method.

public: // interface

  /// Create a new EnzoMethodGravityBiCGStab object
  EnzoMethodGravityBiCGStab(const FieldDescr* field_descr,
			    int rank,
			    double grav_const,
			    int iter_max, 
			    double res_tol,
			    int monitor_iter,
			    bool is_singular,
			    bool diag_precon);

  EnzoMethodGravityBiCGStab() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityBiCGStab);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityBiCGStab(CkMigrateMessage* m) {}

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er& p) {

    // NOTE: change this function whenever attributes change
    TRACEPUP;

    Method::pup(p);

    p | A_;
    p | M_;

    p | is_singular_;
    p | rank_;
    p | grav_const_;
    p | iter_max_;
    p | res_tol_;
    p | monitor_iter_;

    p | idensity_;
    p | ipotential_;
    p | ib_;
    p | ix_;
    p | ir_;
    p | ir0_;
    p | ip_;
    p | iy_;
    p | iv_;
    p | iq_;
    p | iu_;

    p | nx_;
    p | ny_;
    p | nz_;

    p | mx_;
    p | my_;
    p | mz_;

    p | gx_;
    p | gy_;
    p | gz_;

    p | rho0_;
    p | beta_d_;
    p | beta_n_;
    p | vr0_;
    p | omega_d_;
    p | omega_n_;
    p | rr_;

    p | iter_;
    p | beta_;
    p | err_;
    p | err_min_;
    p | err_max_;
    p | alpha_;
    p | omega_;
    p | bs_;
    p | bc_;

    p | id_refresh_matvec_;

  }

  
  /// Enter solver for the gravitational potential
  virtual void compute(Block* block) throw();

  /// Name for the solver
  virtual std::string name() throw() { return "gravity_bicgstab"; }

  /// Shifts RHS and sets initial vectors R, R0, and P
  template<class T> void gravity_bicgstab_start_2(EnzoBlock* enzo_block) throw();

  /// Entry into BiCGStab iteration loop
  template<class T> void gravity_bicgstab_loop_0(EnzoBlock* enzo_block) throw();

  /// First preconditioner solve, begins refresh on Y
  template<class T> void gravity_bicgstab_loop_2(EnzoBlock* enzo_block) throw();

  /// First matrix-vector product, begins dot-product DOT(V,R0)
  template<class T> void gravity_bicgstab_loop_4(EnzoBlock* enzo_block) throw();

  /// First vector updates, begins refresh on Q
  template<class T> void gravity_bicgstab_loop_6(EnzoBlock* enzo_block) throw();

  /// Second preconditioner solve, begins refresh on Y
  template<class T> void gravity_bicgstab_loop_8(EnzoBlock* enzo_block) throw();

  /// Second matrix-vector product, begins dot-products DOT(U,U) and DOT(U,Q)
  template<class T> void gravity_bicgstab_loop_10(EnzoBlock* enzo_block) throw();

  /// Second vector updates, begins dot-products DOT(R,R) and DOT(R,R0)
  template<class T> void gravity_bicgstab_loop_12(EnzoBlock* enzo_block) throw();

  /// Updates search direction, begins update on iteration counter
  template<class T> void gravity_bicgstab_loop_14(EnzoBlock* enzo_block) throw();

  /// End of iteration
  template<class T> void gravity_bicgstab_end(EnzoBlock* enzo_block, int retval) throw();

  /// Set bs_ (B sum) by EnzoBlock after reduction
  void set_bs(long double bs) throw() { bs_ = bs; }

  /// Set bc_ (B count) by EnzoBlock after reduction
  void set_bc(long double bc) throw() { bc_ = bc; }

  /// Set rho0_ by EnzoBlock after reduction
  void set_rho0(long double rho0) throw() { rho0_ = rho0; }

  /// Set beta_d_ by EnzoBlock after reduction
  void set_beta_d(long double beta_d) throw() { beta_d_ = beta_d; }

  /// Set vr0_ by EnzoBlock after reduction
  void set_vr0(long double vr0) throw() { vr0_ = vr0; }

  /// Set omega_d_ by EnzoBlock after reduction
  void set_omega_d(long double omega_d) throw() { omega_d_ = omega_d; }

  /// Set omega_n_ by EnzoBlock after reduction
  void set_omega_n(long double omega_n) throw() { omega_n_ = omega_n; }

  /// Set rr_ by EnzoBlock after reduction
  void set_rr(long double rr) throw() { rr_ = rr; }

  /// Set beta_n_ by EnzoBlock after reduction
  void set_beta_n(long double beta_n) throw() { beta_n_ = beta_n; }

  /// Set iter_ by EnzoBlock after reduction
  void set_iter(int iter) throw() { iter_ = iter; }

protected: // methods

  void monitor_output_(EnzoBlock * enzo_block) throw();

  template<class T> void compute_(EnzoBlock * enzo_block) throw();

  /// Compute local contribution to inner-product X*Y
  template<class T> long double dot_(const T* X, const T* Y) const throw();

  template<class T> void zaxpy_(T* Z, double a, const T* X, const T* Y) const throw();
  
  /// Compute local sum of vector elements X_i
  template<class T> long double sum_(const T* X) const throw();

  /// return the number of elements of the vector X
  long double count_() const throw();
  
protected: // attributes

  /// Matrix
  Matrix* A_;

  /// Preconditioner
  Matrix* M_;

  /// Whether you need to subtract of the nullspace of A from b, e.g. fully
  /// periodic or Neumann problems
  bool is_singular_;

  /// Dimensionality of the problem
  int rank_;

  /// Gas constant, e.g. 6.67384e-8 [cgs]
  double grav_const_;

  /// Maximum number of BiCGStab iterations
  int iter_max_;

  /// Convergence tolerance on the residual reduction rz_ / rz0_
  double res_tol_;

  /// How often to display progress
  int monitor_iter_;

  /// Initial residual
  long double rho0_;

  /// Current error
  long double err_;

  /// Minimum error
  long double err_min_;

  /// Maximum error
  long double err_max_;

  /// Density and potential field id's
  int idensity_;
  int ipotential_;

  /// BiCGStab vector id's
  int ib_;
  int ix_;
  int ir_;
  int ir0_;
  int ip_;
  int iy_;
  int iv_;
  int iq_;
  int iu_;

  /// Block field attributes
  int nx_, ny_, nz_;
  int mx_, my_, mz_;
  int gx_, gy_, gz_;

  /// Current BiCGStab iteration
  int iter_;

  /// dot(R_i,R_0)
  long double beta_d_;

  /// dot(R_i+1,R_0)
  long double beta_n_;

  /// dot(V_i,R_0)
  long double vr0_;

  /// dot(U_i,U_i)
  long double omega_d_;

  /// dot(U_i,Q_i)
  long double omega_n_;

  /// dot(R_i,R_i)
  long double rr_;

  /// beta_n_ / beta_d_
  long double beta_;

  /// omega_n_ / omega_d_
  long double omega_;

  /// beta_d_ / vr0_
  long double alpha_;

  /// sum of elements B(i) for singular systems
  long double bs_;

  /// count of elements B(i) for singular systems
  long double bc_;

  /// matvec refresh index
  int id_refresh_matvec_;

};

#endif /* ENZO_ENZO_METHOD_GRAVITY_BICGSTAB_HPP */
