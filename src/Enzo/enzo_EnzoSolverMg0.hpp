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
  (std::string name,
   int monitor_iter,
   int restart_cycle,
   int rank,
   int iter_max,
   double res_tol,
   int index_smooth_pre,
   int index_solve_coarse,
   int index_smooth_post,
   int index_smooth_last,
   Restrict * restrict,
   Prolong * prolong,
   int min_level,
   int max_level,
   int min_level_coarse,
   int max_level_coarse,
   bool is_unigrid);

  EnzoSolverMg0() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverMg0);
  
  /// Charm++ PUP::able migration constructor
  EnzoSolverMg0 (CkMigrateMessage *m)
    :  Solver(m),
       A_(NULL),
       index_smooth_pre_(-1),
       index_solve_coarse_(-1),
       index_smooth_post_(-1),
       index_smooth_last_(-1),
       restrict_(NULL),
       prolong_(NULL),
       rank_(0),
       iter_max_(0),
       res_tol_(0),
       i_sync_restrict_(-1),i_sync_prolong_(-1),i_iter_(-1),i_msg_(-1),
       ib_(0), ic_(0), ir_(0), ix_(0),
       mx_(0),my_(0),mz_(0),
       gx_(0),gy_(0),gz_(0),
       bs_(0), bc_(0), rr_(0), rr_local_(0), rr0_(0)
  {}

  /// Destructor
  virtual ~EnzoSolverMg0 () throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    // NOTE: change this function whenever attributes change

    TRACEPUP;

    Solver::pup(p);

    //    p | A_;
    p | index_smooth_pre_;
    p | index_solve_coarse_;
    p | index_smooth_post_;
    p | index_smooth_last_;
    p | restrict_;
    p | prolong_;
    p | rank_;
    p | iter_max_;
    p | res_tol_;

    p | i_sync_restrict_;
    p | i_sync_prolong_;
    p | i_iter_;
    p | i_msg_;
    
    p | ib_;
    p | ic_;
    p | ir_;
    p | ix_;

    p | mx_;
    p | my_;
    p | mz_;
    p | gx_;
    p | gy_;
    p | gz_;

    p | bc_;
    p | bs_;

    p | rr_;
    p | rr_local_;
    p | rr0_;
  }

  /// Solve the linear system 
  virtual void apply ( std::shared_ptr<Matrix> A, int ix, int ib,
		       Block * block) throw();

  /// Type of this solver
  virtual std::string type() const { return "mg0"; }

  void compute_shift_(EnzoBlock * enzo_block,long double * reduce) throw();
  void compute_correction(EnzoBlock * enzo_block) throw();

  /// Compute the residual
  void compute_residual_(EnzoBlock *) throw();

  /// Pack and unpack residual for restricting to parent
  FieldMsg * pack_residual_ (EnzoBlock *) throw();
  void unpack_residual_(EnzoBlock *, FieldMsg *) throw();

  /// Pack and unpack correction for prolonging to child
  FieldMsg * pack_correction_(EnzoBlock * enzo_block, int ic3[3]) throw();
  void unpack_correction_(EnzoBlock *, FieldMsg *) throw();
  
  /// Restrict residual to coarser Block
  void restrict(EnzoBlock * enzo_block) throw();

  /// Restrict residual to parent
  void restrict_send(EnzoBlock * enzo_block) throw();
  void restrict_recv(EnzoBlock * enzo_block,
		     FieldMsg * field_message) throw();

  // /// Access the Restrict operator by EnzoBlock
  // Restrict * restrict() { return restrict_; }

  // /// Access the Prolong operator by EnzoBlock
  // Prolong * prolong() { return prolong_; }

  /// Call coarse solver--must be called by all blocks
  void call_coarse_solver(EnzoBlock * enzo_block) throw();
  /// Call pre-smoother--must be called by all blocks (or not at all)
  void call_pre_smoother(EnzoBlock * enzo_block) throw();
  /// Call post-smoother--must be called by all blocks (or not at all)
  void call_post_smoother(EnzoBlock * enzo_block) throw();

  /// Begin the prolongation phase
  void prolong(EnzoBlock * enzo_block) throw();

  /// Prolong the correction C to the next-finer level
  void prolong_recv(EnzoBlock * enzo_block,
		    FieldMsg * field_message) throw();

  /// Apply post-smoothing to the current level
  void post_smooth(EnzoBlock * enzo_block) throw();

  void set_bs(long double bs) throw() { bs_ = bs; }
  void set_bc(long double bc) throw() { bc_ = bc; }
  void set_rr_local(long double rr) throw() { rr_local_ = rr; }
  void set_rr(long double rr) throw() { rr_ = rr; }
  void set_rr0(long double rr0) throw() { rr0_ = rr0; }

  long double rr_local() throw() { return rr_local_; }
  long double rr() throw() { return rr_; }

  void begin_solve(EnzoBlock * enzo_block,
		   CkReductionMsg *msg) throw();

  void end_cycle(EnzoBlock * enzo_block) throw();
  
  void print()
  {
    //    CkPrintf (" A_ = %p\n",A_);
    CkPrintf (" index_smooth_pre_ = %d\n",index_smooth_pre_);
    CkPrintf (" index_solve_coarse_ = %d\n",index_solve_coarse_);
    CkPrintf (" index_smooth_post_ = %d\n",index_smooth_post_);
    CkPrintf (" index_smooth_last_ = %d\n",index_smooth_last_);
    CkPrintf (" restrict_ = %p\n",restrict_);
    CkPrintf (" prolong_ = %p\n",prolong_);
    CkPrintf (" rank_ = %d\n",rank_);
    CkPrintf (" iter_max_ = %d\n",iter_max_);
    CkPrintf (" res_tol_ = %g\n",res_tol_);
    CkPrintf (" i_sync_restrict_ = %g\n",i_sync_restrict_);
    CkPrintf (" i_sync_prolong_ = %g\n",i_sync_prolong_);
    CkPrintf (" i_iter_ = %g\n",i_iter_);
    CkPrintf (" i_msg_ = %g\n",i_msg_);
    CkPrintf (" ib_ = %d\n",ib_);
    CkPrintf (" ic_ = %d\n",ic_);
    CkPrintf (" ir_ = %d\n",ir_);
    CkPrintf (" ix_ = %d\n",ix_);
    CkPrintf (" mx_,my_,mz_ = %d %d %d\n",mx_,my_,mz_);
    CkPrintf (" gx_,gy_,gz_ = %d %d %d\n",gx_,gy_,gz_);
    CkPrintf (" bs_ = %g\n",bs_);
    CkPrintf (" bc_ = %g\n",bc_);
    CkPrintf (" rr_ = %g\n",rr_);
    CkPrintf (" rr_local_ = %g\n",rr_local_);
    CkPrintf (" rr0_ = %g\n",rr0_);
  }

  /// Exit the solver
  void end(Block * block);


  /// Access the prolong Sync Scalar value for the Block
  Sync * psync_prolong(Block * block)
  {
    ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
    ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
    return scalar_data->value(scalar_descr,i_sync_prolong_);
  }

  /// Access the restrict Sync Scalar value for the Block
  Sync * psync_restrict(Block * block)
  {
    ScalarData<Sync> * scalar_data = block->data()->scalar_data_sync();
    ScalarDescr *      scalar_descr = cello::scalar_descr_sync();
    return scalar_data->value(scalar_descr,i_sync_restrict_);
  }

  /// Access the iter Scalar value for the Block
  int * piter(Block * block)
  {
    ScalarData<int> * scalar_data = block->data()->scalar_data_int();
    ScalarDescr *      scalar_descr = cello::scalar_descr_int();
    return scalar_data->value(scalar_descr,i_iter_);
  }
  
  /// Access the msg Scalar value for the Block
  FieldMsg ** pmsg(Block * block)
  {
    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();
    ScalarDescr *        scalar_descr = cello::scalar_descr_void();
    return (FieldMsg **)scalar_data->value(scalar_descr,i_msg_);
  }
  
protected: // methods

  void enter_solver_(EnzoBlock * enzo_block) throw();

  void begin_cycle_(EnzoBlock * enzo_block) throw();
  
  /// Prolong the correction C to the next-finer level
  void prolong_send_(EnzoBlock * enzo_block) throw();

  bool is_converged_(EnzoBlock * enzo_block) const;
  bool is_diverged_(EnzoBlock * enzo_block) const;

  /// Shift RHS if needed for singular problems
  void do_shift_b_(EnzoBlock *, CkReductionMsg *) throw();
  
  /// Allocate temporary Fields
  void allocate_temporary_(Field field, Block * block = NULL)
  {
    field.allocate_temporary(ir_);
    field.allocate_temporary(ic_);
  }
	      
  /// Dellocate temporary Fields

  void deallocate_temporary_(Field field, Block * block = NULL)
  {
    field.deallocate_temporary(ir_);
    field.deallocate_temporary(ic_);
  }
  
protected: // attributes

  /// Matrix
  std::shared_ptr<Matrix> A_;

  /// Solver for pre-smoother
  int index_smooth_pre_;

  /// Solver for coarse solver
  int index_solve_coarse_;

  /// Solver for post smoother
  int index_smooth_post_;

  /// Solver for final smoothing of solution
  int index_smooth_last_;

  /// Restriction
  Restrict * restrict_;

  /// Prolongation
  Prolong * prolong_;

  /// Dimensionality of the problem
  int rank_;

  /// Maximum number of MG iterations
  int iter_max_;

  /// Convergence tolerance on the residual reduction rr_ / rr0_
  double res_tol_;

  /// MG scalar id's
  int i_sync_restrict_;
  int i_sync_prolong_;
  int i_iter_;
  int i_msg_;

  /// MG vector id's
  int ib_;
  int ic_;
  int ir_;
  int ix_;

  /// Block field attributes
  int mx_,my_,mz_;
  int gx_,gy_,gz_;

  /// scalars used for projections of singular systems
  long double bs_;
  long double bc_;

  /// Current and initial residual norm R'*R
  long double rr_;
  long double rr_local_;
  long double rr0_;

  /// Min and Max level for coarse-grid solver (in case coarse grid solver is multilevel)
  int min_level_coarse_;
  int max_level_coarse_;
};

#endif /* ENZO_ENZO_SOLVER_GRAVITY_MG0_HPP */
