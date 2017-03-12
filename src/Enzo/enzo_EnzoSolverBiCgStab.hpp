// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverBiCgStab.hpp
/// @author   Daniel R. Reynolds (reynolds@smu.edu)
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2014-10-21 17:25:40
/// @brief    [\ref Enzo] Declaration of EnzoSolverBiCgStab
///
/// Bicongugate gradient stabilized solver (BiCgStab) for solving 
/// linear systems on field data.

#ifndef ENZO_ENZO_SOLVER_BICGSTAB_HPP
#define ENZO_ENZO_SOLVER_BICGSTAB_HPP

class EnzoSolverBiCgStab : public Solver {

  /// @class    EnzoSolverBiCgStab
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] This class implements the BiCgStab Krylov
  /// linear solver.  Alone, this is more applicable to smaller
  /// problems since the solver doesn't scale as well as some other
  /// solvers (FFT, MG, etc.) for larger problems.  Alternately, a
  /// more scalable solver may be combined as a preconditioner for a
  /// robust and scalable overall solver.

public: // interface

  /// normal constructor
  EnzoSolverBiCgStab(const FieldDescr* field_descr,
		     int monitor_iter,
		     int rank,
		     int iter_max, 
		     double res_tol,
		     int min_level,
		     int max_level,
		     int index_precon);

  /// default constructor
  EnzoSolverBiCgStab()
    : Solver(),
      A_(NULL),
      index_precon_(-1),
      first_call_(true),
      rank_(0),
      iter_max_(0), 
      res_tol_(0.0),
      rho0_(0), err_(0), err0_(0), err_min_(0), err_max_(0),
      idensity_(0),  ipotential_(0),
      ib_(0), ix_(0), ir_(0), ir0_(0), ip_(0), 
      iy_(0), iv_(0), iq_(0), iu_(0),
      nx_(0), ny_(0), nz_(0),
      mx_(0), my_(0), mz_(0),
      gx_(0), gy_(0), gz_(0),
      iter_(0),
      beta_d_(0), beta_n_(0), beta_(0), 
      omega_d_(0), omega_n_(0), omega_(0), 
      vr0_(0), rr_(0), alpha_(0),
      bs_(0.0),bc_(0.0),
      ys_(0.0),vs_(0.0),us_(0.0)
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverBiCgStab);
  
  /// destructor
  ~EnzoSolverBiCgStab() throw();
  
  /// Charm++ PUP::able migration constructor
  EnzoSolverBiCgStab(CkMigrateMessage* m)
    : Solver(m),
      A_(NULL),
      first_call_(true),
      rank_(0),
      iter_max_(0), 
      res_tol_(0.0),
      rho0_(0), err_(0), err0_(0),err_min_(0), err_max_(0),
      idensity_(0),  ipotential_(0),
      ib_(0), ix_(0), ir_(0), ir0_(0), ip_(0), 
      iy_(0), iv_(0), iq_(0), iu_(0),
      nx_(0), ny_(0), nz_(0),
      mx_(0), my_(0), mz_(0),
      gx_(0), gy_(0), gz_(0),
      iter_(0),
      beta_d_(0), beta_n_(0), beta_(0), 
      omega_d_(0), omega_n_(0), omega_(0), 
      vr0_(0), rr_(0), alpha_(0),
      bs_(0.0),bc_(0.0),
      ys_(0.0),vs_(0.0),us_(0.0)
  {}

  /// Charm++ Pack / Unpack function
  void pup(PUP::er& p) {

    // JB NOTE: change this function whenever attributes change
    TRACEPUP;

    Solver::pup(p);

    p | A_;
    p | index_precon_;
    
    p | first_call_;
    p | rank_;
    p | iter_max_;
    p | res_tol_;

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
    p | err0_;
    p | err_min_;
    p | err_max_;
    p | alpha_;
    p | omega_;
    p | bs_;
    p | bc_;
    p | ys_;
    p | vs_;
    p | us_;

  }

  
  /// Main solver entry routine
  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw();

  /// Name to call the solver within Enzo-P
  virtual std::string name() const
  { return "bicgstab"; }

  /// Projects RHS and sets initial vectors R, R0, and P
  template<class T> void start_2(EnzoBlock* enzo_block) throw();

  /// Entry into BiCgStab iteration loop, begins refresh on P
  template<class T> void loop_0(EnzoBlock* enzo_block) throw();

  /// First preconditioner solve
  template<class T> void loop_2(EnzoBlock* enzo_block) throw();

  /// Return from preconditioner solve, begins refresh on Y
  template<class T> void loop_25(EnzoBlock* enzo_block) throw();

  /// First matrix-vector product, begins DOT(V,R0) and projection of
  /// Y and V
  template<class T> void loop_4(EnzoBlock* enzo_block) throw();

  /// Shifts Y and V, begins, first vector updates, begins refresh on Q
  template<class T> void loop_6(EnzoBlock* enzo_block) throw();

  /// Second preconditioner solve, begins refresh on Y
  template<class T> void loop_8(EnzoBlock* enzo_block) throw();

  /// Return from preconditioner solve, begins refresh on Y
  template<class T> void loop_85(EnzoBlock* enzo_block) throw();

  /// Second matrix-vector product, begins DOT(U,U), DOT(U,Q) and
  /// projection of Y and U
  template<class T> void loop_10(EnzoBlock* enzo_block) throw();

  /// Shifts Y and U, second vector updates, begins DOT(R,R) and
  /// DOT(R,R0)
  template<class T> void loop_12(EnzoBlock* enzo_block) throw();

  /// Updates search direction, begins update on iteration counter
  template<class T> void loop_14(EnzoBlock* enzo_block) throw();

  /// End of iteration
  template<class T> void end(EnzoBlock* enzo_block, int retval) throw();

  /// Compute accelerations
  template<class T> void acc(EnzoBlock* enzo_block) throw();

  /// Exit the solver
  template<class T> void exit(EnzoBlock* enzo_block) throw();

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

  /// internal routine to handle actual start to solver
  template<class T> void compute_(EnzoBlock * enzo_block) throw();

  /// Compute local contribution to inner-product X*Y
  template<class T> long double dot_(const T* X, const T* Y) const throw();

  /// Perform local vector update: Z = A*X + Y
  template<class T> void zaxpy_(T* Z, double a, const T* X, const T* Y) const throw();
  
  /// Compute local sum of vector elements X_i
  template<class T> long double sum_(const T* X) const throw();

protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Matrix
  Matrix* A_;

  /// Preconditioner (-1 if none)
  int index_precon_;

  /// Whether this is the first call to the solver
  bool first_call_;
  
  /// Dimensionality of the problem
  int rank_;

  /// Maximum number of allowed BiCgStab iterations
  int iter_max_;

  /// Convergence tolerance on the relative residual
  double res_tol_;

  /// Initial residual
  long double rho0_;

  /// Current error
  long double err_;

  /// Initial error
  long double err0_;

  /// Minimum error (all iterations so far)
  long double err_min_;

  /// Maximum error (all iterations so far)
  long double err_max_;

  /// Density and potential field id's (for RHS and solution)
  int idensity_;
  int ipotential_;

  /// BiCgStab vector id's
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

  /// Current BiCgStab iteration
  int iter_;

  /// scalars used within BiCgStab iteration
  long double beta_d_;
  long double beta_n_;
  long double beta_;
  long double omega_d_;
  long double omega_n_;
  long double omega_;
  long double vr0_;
  long double rr_;
  long double alpha_;

  /// scalars used for projections of singular systems
  long double bs_;
  long double bc_;
  long double ys_;
  long double vs_;
  long double us_;

};

#endif /* ENZO_ENZO_SOLVER_BICGSTAB_HPP */
