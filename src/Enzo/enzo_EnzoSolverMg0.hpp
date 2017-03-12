// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverMg0.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2015-06-02
/// @brief    [\ref Enzo] Declaration of EnzoSolverMg0
///
/// Multigrid for solving a linear system on the root-level grid only

#ifndef ENZO_ENZO_SOLVER_MG0_HPP
#define ENZO_ENZO_SOLVER_MG0_HPP

class EnzoSolverMg0 : public Solver {

  /// @class    EnzoSolverMg0
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Multigrid on the root-level grid.  For use either
  /// as a Gravity solver on non-adaptive problems, or as a preconditioner
  /// for a Krylov subspace solver, as in Dan Reynold's HG solver.

public: // interface

  /// Create a new EnzoSolverMg0 object
  EnzoSolverMg0
  (const FieldDescr * field_descr,
   int monitor_iter,
   int rank,
   int iter_max,
   int index_smooth_pre,
   int index_solve_coarse,
   int index_smooth_post,
   Restrict * restrict,
   Prolong * prolong,
   int min_level,
   int max_level);

  EnzoSolverMg0() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverMg0);
  
  /// Charm++ PUP::able migration constructor
  EnzoSolverMg0 (CkMigrateMessage *m)
    :  A_(NULL),
       index_smooth_pre_(-1),
       index_solve_coarse_(-1),
       index_smooth_post_(-1),
       restrict_(NULL),
       prolong_(NULL),
       rank_(0),
       iter_max_(0), 
       monitor_iter_(0),
       ib_(0), ic_(0), ir_(0), ix_(0),
       mx_(0),my_(0),mz_(0),
       nx_(0),ny_(0),nz_(0),
       gx_(0),gy_(0),gz_(0),
       bs_(0), bc_(0)
  {}

  /// Destructor
  ~EnzoSolverMg0 () throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    // NOTE: change this function whenever attributes change

    TRACEPUP;

    Solver::pup(p);

    p | A_;
    p | index_smooth_pre_;
    p | index_solve_coarse_;
    p | index_smooth_post_;
    p | restrict_;
    p | prolong_;
    p | rank_;
    p | iter_max_;
    p | monitor_iter_;

    p | ib_;
    p | ic_;
    p | ir_;
    p | ix_;

    p | mx_;
    p | my_;
    p | mz_;
    p | nx_;
    p | ny_;
    p | nz_;
    p | gx_;
    p | gy_;
    p | gz_;

    p | bc_;
    p | bs_;
  }

  /// Solve the linear system 
  virtual void apply ( Matrix * A, int ix, int ib, Block * block) throw();

  virtual std::string name () const
  { return "mg0"; }

  void compute_correction(EnzoBlock * enzo_block) throw();

  /// Apply pre-smoothing on the current level
  template <class T>
  void pre_smooth(EnzoBlock * enzo_block) throw();

  /// Restrict residual to parent
  template <class T>
  void restrict_send(EnzoBlock * enzo_block) throw();
  template <class T>
  void restrict_recv(EnzoBlock * enzo_block,
		      FieldMsg * field_message) throw();

  /// Access the Restrict operator by EnzoBlock
  Restrict * restrict() { return restrict_; }

  /// Access the Prolong operator by EnzoBlock
  Prolong * prolong() { return prolong_; }
  
  /// Solve the coarse-grid equation A*C = R
  template <class T>
  void solve_coarse(EnzoBlock * enzo_block) throw();

  /// Prolong the correction C to the next-finer level
  template <class T>
  void prolong_recv(EnzoBlock * enzo_block,
		    FieldMsg * field_message) throw();

  /// Apply post-smoothing to the current level
  template <class T>
  void post_smooth(EnzoBlock * enzo_block) throw();

  void set_bs(long double bs) throw() { bs_ = bs; }
  void set_bc(long double bc) throw() { bc_ = bc; }

  template <class T>
  void begin_solve(EnzoBlock * enzo_block) throw();

  template <class T>
  void end_cycle(EnzoBlock * enzo_block) throw();
  
  void print()
  {
    CkPrintf (" A_ = %p\n",A_);
    CkPrintf (" index_smooth_pre_ = %d\n",index_smooth_pre_);
    CkPrintf (" index_solve_coarse_ = %d\n",index_solve_coarse_);
    CkPrintf (" index_smooth_post_ = %d\n",index_smooth_post_);
    CkPrintf (" restrict_ = %p\n",restrict_);
    CkPrintf (" prolong_ = %p\n",prolong_);
    CkPrintf (" rank_ = %d\n",rank_);
    CkPrintf (" iter_max_ = %d\n",iter_max_);
    CkPrintf (" monitor_iter_ = %d\n",monitor_iter_);
    CkPrintf (" rr0_ = %g\n",rr0_);
    CkPrintf (" ib_ = %d\n",ib_);
    CkPrintf (" ic_ = %d\n",ic_);
    CkPrintf (" ir_ = %d\n",ir_);
    CkPrintf (" ix_ = %d\n",ix_);
    CkPrintf (" mx_,my_,mz_ = %d %d %d\n",mx_,my_,mz_);
    CkPrintf (" nx_,ny_,nz_ = %d %d %d\n",nx_,ny_,nz_);
    CkPrintf (" gx_,gy_,gz_ = %d %d %d\n",gx_,gy_,gz_);
    CkPrintf (" bs_ = %g\n",bs_);
    CkPrintf (" bc_ = %g\n",bc_);
  }
protected: // methods

  template <class T>
  void enter_solver_(EnzoBlock * enzo_block) throw();

  template <class T>
  void begin_cycle_(EnzoBlock * enzo_block) throw();
  
  /// Prolong the correction C to the next-finer level
  template <class T>
  void prolong_send_(EnzoBlock * enzo_block) throw();

  void monitor_output_(EnzoBlock * enzo_block) throw();

  bool is_converged_(EnzoBlock * enzo_block) const;

protected: // attributes

  /// Matrix
  Matrix * A_;

  /// Solver for pre-smoother
  int index_smooth_pre_;

  /// Solver for coarse solver
  int index_solve_coarse_;

  /// Solver for post smoother
  int index_smooth_post_;

  /// Restriction
  Restrict * restrict_;

  /// Prolongation
  Prolong * prolong_;

  /// Dimensionality of the problem
  int rank_;

  /// Maximum number of MG iterations
  int iter_max_;

  /// How often to display progress
  int monitor_iter_;

  /// Initial residual
  long double rr0_;

  /// MG vector id's
  int ib_;
  int ic_;
  int ir_;
  int ix_;

  /// Block field attributes
  int mx_,my_,mz_;
  int nx_,ny_,nz_;
  int gx_,gy_,gz_;

  /// scalars used for projections of singular systems
  long double bs_;
  long double bc_;
};

#endif /* ENZO_ENZO_SOLVER_GRAVITY_MG0_HPP */
