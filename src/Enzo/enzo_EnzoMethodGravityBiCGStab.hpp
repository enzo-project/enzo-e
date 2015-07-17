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

  /// default constructor
  EnzoMethodGravityBiCGStab() {};

  /// normal constructor
  EnzoMethodGravityBiCGStab(FieldDescr* field_descr,
			    int rank,
			    double grav_const,
			    int iter_max, 
			    double res_tol,
			    int monitor_iter,
			    bool is_singular,
			    bool diag_precon);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityBiCGStab);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityBiCGStab(CkMigrateMessage* m) {}

  /// Charm++ Pack / Unpack function
  void pup(PUP::er& p) {

    // JB NOTE: change this function whenever attributes change
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
    p | ys_;
    p | vs_;
    p | us_;

    p | id_refresh_P_;
    p | id_refresh_Q_;
    p | id_refresh_Y_;

  }

  
  /// Main solver entry routine
  virtual void compute(Block* block) throw();

  /// Name to call the solver within Enzo-P
  virtual std::string name() throw() { return "gravity_bicgstab"; }

  /// Projects RHS and sets initial vectors R, R0, and P
  template<class T> void start_2(EnzoBlock* enzo_block) throw();

  /// Entry into BiCGStab iteration loop, begins refresh on P
  template<class T> void loop_0(EnzoBlock* enzo_block) throw();

  /// First preconditioner solve, begins refresh on Y
  template<class T> void loop_2(EnzoBlock* enzo_block) throw();

  /// First matrix-vector product, begins DOT(V,R0) and projection of Y and V
  template<class T> void loop_4(EnzoBlock* enzo_block) throw();

  /// Shifts Y and V, begins, first vector updates, begins refresh on Q
  template<class T> void loop_6(EnzoBlock* enzo_block) throw();

  /// Second preconditioner solve, begins refresh on Y
  template<class T> void loop_8(EnzoBlock* enzo_block) throw();

  /// Second matrix-vector product, begins DOT(U,U), DOT(U,Q) and projection of Y and U
  template<class T> void loop_10(EnzoBlock* enzo_block) throw();

  /// Shifts Y and U, second vector updates, begins DOT(R,R) and DOT(R,R0)
  template<class T> void loop_12(EnzoBlock* enzo_block) throw();

  /// Updates search direction, begins update on iteration counter
  template<class T> void loop_14(EnzoBlock* enzo_block) throw();

  /// End of iteration
  template<class T> void end(EnzoBlock* enzo_block, int retval) throw();

  /// Set routines for use by EnzoBlock after reductions
  void set_bs(long double bs) throw() { bs_ = bs; }
  void set_bc(long double bc) throw() { bc_ = bc; }
  void set_ys(long double ys) throw() { ys_ = ys; }
  void set_vs(long double vs) throw() { vs_ = vs; }
  void set_us(long double us) throw() { us_ = us; }
  void set_rho0(long double rho0) throw() { rho0_ = rho0; }
  void set_beta_d(long double beta_d) throw() { beta_d_ = beta_d; }
  void set_vr0(long double vr0) throw() { vr0_ = vr0; }
  void set_omega_d(long double omega_d) throw() { omega_d_ = omega_d; }
  void set_omega_n(long double omega_n) throw() { omega_n_ = omega_n; }
  void set_rr(long double rr) throw() { rr_ = rr; }
  void set_beta_n(long double beta_n) throw() { beta_n_ = beta_n; }
  void set_iter(int iter) throw() { iter_ = iter; }

protected: // methods

  /// internal routine to report solver progress to stdout
  void monitor_output_(EnzoBlock * enzo_block) throw();

  /// internal routine to handle actual start to solver
  template<class T> void compute_(EnzoBlock * enzo_block) throw();

  /// Compute local contribution to inner-product X*Y
  template<class T> long double dot_(const T* X, const T* Y) const throw();

  /// Perform local vector update: Z = A*X + Y
  template<class T> void zaxpy_(T* Z, double a, const T* X, const T* Y) const throw();
  
  /// Compute local sum of vector elements X_i
  template<class T> long double sum_(const T* X) const throw();

  /// return the number of elements of the vector on this block
  long double count_() const throw();
  
protected: // attributes

  /// Matrix
  Matrix* A_;

  /// Preconditioner
  Matrix* M_;

  /// Whether you need to project b into R(A), e.g. fully periodic or Neumann problems
  bool is_singular_;

  /// Dimensionality of the problem
  int rank_;

  /// Gravitational constant, e.g. 6.67384e-8 [cgs]
  double grav_const_;

  /// Maximum number of allowed BiCGStab iterations
  int iter_max_;

  /// Convergence tolerance on the relative residual
  double res_tol_;

  /// How often to display progress to stdout
  int monitor_iter_;

  /// Initial residual
  long double rho0_;

  /// Current error
  long double err_;

  /// Minimum error (all iterations so far)
  long double err_min_;

  /// Maximum error (all iterations so far)
  long double err_max_;

  /// Density and potential field id's (for RHS and solution)
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
  int nx_, ny_, nz_;   /// active block size
  int mx_, my_, mz_;   /// total block size
  int gx_, gy_, gz_;   /// ghost zones

  /// Current BiCGStab iteration
  int iter_;

  /// scalars used within BiCGStab iteration
  long double beta_d_;
  long double beta_n_;
  long double beta_;
  long double omega_d_;
  long double omega_n_;
  long double omega_;
  long double vr0_;
  long double rr_;
  long double alpha_;

  /// scalars used for projections of singular gravity systems
  long double bs_;
  long double bc_;
  long double ys_;
  long double vs_;
  long double us_;

  /// refresh indices
  int id_refresh_P_;
  int id_refresh_Q_;
  int id_refresh_Y_;

};

#endif /* ENZO_ENZO_METHOD_GRAVITY_BICGSTAB_HPP */
