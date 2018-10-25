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
  EnzoSolverBiCgStab(std::string name,
		     std::string field_x,
		     std::string field_b,
		     int monitor_iter,
		     int reuse_solution,
		     int rank,
		     int iter_max, 
		     double res_tol,
		     int min_level,
		     int max_level,
		     int index_precon,
		     int solve_type);

  /// default constructor
  EnzoSolverBiCgStab()
    : Solver(),
      alpha_(0), beta_n_(0), beta_d_(0),   omega_(0),
      rr_(0), r0s_(0.0), c_(0.0), bnorm_(0.0),
      rho0_(0), err_(0), err0_(0), err_min_(0), err_max_(0),
      res_tol_(0.0),
      A_(NULL),
      index_precon_(-1),
      iter_max_(0), 
      ir_(-1), ir0_(-1), ip_(-1), 
      iy_(-1), iv_(-1), iq_(-1), iu_(-1),
      nx_(0), ny_(0), nz_(0),
      m_(0), mx_(0), my_(0), mz_(0),
      gx_(0), gy_(0), gz_(0),
      iter_(0)
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverBiCgStab);
  
  /// Charm++ PUP::able migration constructor
  EnzoSolverBiCgStab(CkMigrateMessage* m)
    : Solver(m),
      alpha_(0), beta_n_(0), beta_d_(0),   omega_(0),
      rr_(0), r0s_(0.0), c_(0.0), bnorm_(0.0),
      rho0_(0), err_(0), err0_(0), err_min_(0), err_max_(0),
      res_tol_(0.0),
      A_(NULL),
      index_precon_(-1),
      iter_max_(0), 
      ir_(-1), ir0_(-1), ip_(-1), 
      iy_(-1), iv_(-1), iq_(-1), iu_(-1),
      nx_(0), ny_(0), nz_(0),
      m_(0), mx_(0), my_(0), mz_(0),
      gx_(0), gy_(0), gz_(0),
      iter_(0)
  {}

  /// Charm++ Pack / Unpack function
  void pup(PUP::er& p) {

    // JB NOTE: change this function whenever attributes change
    TRACEPUP;

    Solver::pup(p);

    //    p | A_;
    p | index_precon_;
    
    p | iter_max_;
    p | res_tol_;

    p | rho0_;
    p | err_;
    p | err0_;
    p | err_min_;
    p | err_max_;

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

    p | m_;
    p | mx_;
    p | my_;
    p | mz_;

    p | gx_;
    p | gy_;
    p | gz_;

    p | iter_;

    p | alpha_;
    p | beta_n_;
    p | beta_d_;
    p | omega_;

    p | rr_;
    p | r0s_;
    p | c_;
    p | bnorm_;

  }

  
  /// Main solver entry routine
  virtual void apply (std::shared_ptr<Matrix> A, Block * block) throw();

  /// Type of this solver
  virtual std::string type() const { return "bicgstab"; }

  /// Projects RHS and sets initial vectors R, R0, and P
  void start_2(EnzoBlock* enzo_block,
	       CkReductionMsg * msg) throw();

  /// Entry into BiCgStab iteration loop, begins refresh on P
  void loop_0a(EnzoBlock* enzo_block,
	      CkReductionMsg *) throw();
  void loop_0b(EnzoBlock* enzo_block,
	      CkReductionMsg *) throw();
  void loop_0(EnzoBlock* enzo_block) throw();

  /// First preconditioner solve
  void loop_2(EnzoBlock* enzo_block) throw();

  /// Return from preconditioner solve, begins refresh on Y
  void loop_25(EnzoBlock* enzo_block) throw();

  /// First matrix-vector product, begins DOT(V,R0) and projection of
  /// Y and V
  void loop_4(EnzoBlock* enzo_block) throw();

  /// Shifts Y and V, begins, first vector updates, begins refresh on Q
  void loop_6(EnzoBlock* enzo_block, CkReductionMsg *) throw();

  /// Second preconditioner solve, begins refresh on Y
  void loop_8(EnzoBlock* enzo_block) throw();

  /// Return from preconditioner solve, begins refresh on Y
  void loop_85(EnzoBlock* enzo_block) throw();

  /// Second matrix-vector product, begins DOT(U,U), DOT(U,Q) and
  /// projection of Y and U
  void loop_10(EnzoBlock* enzo_block) throw();

  /// Shifts Y and U, second vector updates, begins DOT(R,R) and
  /// DOT(R,R0)
  void loop_12(EnzoBlock* enzo_block, CkReductionMsg * ) throw();

  /// Updates search direction, begins update on iteration counter
  void loop_14(EnzoBlock* enzo_block, CkReductionMsg * ) throw();

  /// End the solve
  void end(EnzoBlock* enzo_block, int retval) throw();

protected: // methods

  /// internal routine to handle actual start to solver
  void compute_(EnzoBlock * enzo_block) throw();

  /// Allocate temporary Fields
  void allocate_temporary_(Block * block)
  {
    Field field = block->data()->field();
    field.allocate_temporary(ir_);
    field.allocate_temporary(ir0_);
    field.allocate_temporary(ip_);
    field.allocate_temporary(iy_);
    field.allocate_temporary(iv_);
    field.allocate_temporary(iq_);
    field.allocate_temporary(iu_);
  }

  /// Dellocate temporary Fields
  void deallocate_temporary_(Block * block)
  {
    Field field = block->data()->field();
    field.deallocate_temporary(ir_);
    field.deallocate_temporary(ir0_);
    field.deallocate_temporary(ip_);
    field.deallocate_temporary(iy_);
    field.deallocate_temporary(iv_);
    field.deallocate_temporary(iq_);
    field.deallocate_temporary(iu_);
  }
  
protected: // attributes

  // NOTE: change pup() function whenever attributes change

  /// scalars used within BiCgStab iteration

  double alpha_;
  double beta_n_;
  double beta_d_;
  double omega_;
  double rr_;
  double r0s_; // sum (R0[i])
  double c_;  // B.length() ("count")
  double bnorm_; // used when reuse_solution

  /// Initial residual
  double rho0_;

  /// Current error
  double err_;

  /// Initial error
  double err0_;

  /// Minimum error (all iterations so far)
  double err_min_;

  /// Maximum error (all iterations so far)
  double err_max_;

  /// Convergence tolerance on the relative residual
  double res_tol_;

  /// Matrix
  std::shared_ptr<Matrix> A_;

  /// Preconditioner (-1 if none)
  int index_precon_;

  /// Maximum number of allowed BiCgStab iterations
  int iter_max_;

  /// BiCgStab vector id's
  int ir_;
  int ir0_;
  int ip_;
  int iy_;
  int iv_;
  int iq_;
  int iu_;

  /// Block field attributes
  int nx_, ny_, nz_;   /// active block size
  int m_;              /// product mx_*my_*mz_ for convenience
  int mx_, my_, mz_;   /// total block size
  int gx_, gy_, gz_;   /// ghost zones

  /// Current BiCgStab iteration
  int iter_;

};

#endif /* ENZO_ENZO_SOLVER_BICGSTAB_HPP */
