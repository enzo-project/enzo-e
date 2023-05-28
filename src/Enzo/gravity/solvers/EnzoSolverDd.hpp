// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSolverDd.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2018-10-01
/// @brief    [\ref Enzo] Declaration of EnzoSolverDd
///
/// Domain decomposition solver

#ifndef ENZO_ENZO_SOLVER_DD_HPP
#define ENZO_ENZO_SOLVER_DD_HPP

class EnzoSolverDd : public Solver {

  /// @class    EnzoSolverDd
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Multigrid on the root-level grid using Mg0, then
  /// BiCgStab in overlapping subdomains defined by root-level Blocks.
  /// An optional final Jacobi step can be applied to smooth the solution
  /// along subdomain boundaries.

public: // interface

  /// Create a new EnzoSolverDd object
  EnzoSolverDd
  (std::string name,
   std::string field_x,
   std::string field_b,
   int monitor_iter,
   int restart_cycle,
   int solve_type,
   int index_prolong,
   int index_restrict,
   int min_level,
   int max_level,
   int index_solve_coarse,
   int index_solve_domain,
   int index_solve_smooth,
   int coarse_level) ;

  EnzoSolverDd() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoSolverDd);
  
  /// Charm++ PUP::able migration constructor
  EnzoSolverDd (CkMigrateMessage *m)
    :  Solver(m),
       A_(NULL),
       index_solve_coarse_(-1),
       index_solve_domain_(-1),
       index_solve_smooth_(-1),
       index_prolong_(-1),
       index_restrict_(-1),
       i_sync_restrict_(-1),
       i_sync_prolong_(-1),
       i_msg_restrict_(),
       i_msg_prolong_(-1),
       ixc_(-1),
       mx_(0),my_(0),mz_(0),
       gx_(0),gy_(0),gz_(0),
       coarse_level_(0)
  {
    for (int i=0; i<cello::num_children(); i++) i_msg_restrict_[i] = -1;
  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    // NOTE: change this function whenever attributes change

    TRACEPUP;

    Solver::pup(p);

    //    p | A_;
    p | index_solve_coarse_;
    p | index_solve_domain_;
    p | index_solve_smooth_;
    p | index_prolong_;
    p | index_restrict_;

    p | i_sync_restrict_;
    p | i_sync_prolong_;
    PUParray(p,i_msg_restrict_,8);
    p | i_msg_prolong_;

    p | ixc_;

    p | mx_;
    p | my_;
    p | mz_;
    p | gx_;
    p | gy_;
    p | gz_;

    p | coarse_level_;

  }

public:  // virtual methods
  
  /// Solve the linear system 
  virtual void apply ( std::shared_ptr<Matrix> A, Block * block) throw();

  /// Type of this solver
  virtual std::string type() const { return "dd"; }

public: // methods

  /// synchronize entering solver for each level
  void begin_solve(EnzoBlock * enzo_block) throw();

  /// Restrict b to coarser Block
  void restrict(EnzoBlock * enzo_block) throw();
  void restrict_send(EnzoBlock * enzo_block) throw();
  void restrict_recv(EnzoBlock * enzo_block,
		     FieldMsg * field_message) throw();

  /// Call coarse solver--must be called by all blocks
  void call_coarse_solver(EnzoBlock * enzo_block) throw();
  
  /// Begin the prolongation phase
  void prolong(EnzoBlock * enzo_block) throw();
  void prolong_send_(EnzoBlock * enzo_block) throw();
  void prolong_recv(EnzoBlock * enzo_block,
		    FieldMsg * field_message) throw();

  /// Call domain solver
  void call_domain_solver(EnzoBlock * enzo_block) throw();
  void continue_after_domain_solve(EnzoBlock * enzo_block) throw();

  /// Call last smoother
  void call_last_smoother(EnzoBlock * enzo_block) throw();
  void continue_after_last_smooth(EnzoBlock * enzo_block) throw();
  
  /// Allocate temporary Fields
  void allocate_temporary_(Block * block)
  {
    Field field = block->data()->field();
    field.allocate_temporary(ixc_);
  }
	      
  /// Dellocate temporary Fields

  void deallocate_temporary_(Block * block)
  {
    Field field = block->data()->field();
    field.deallocate_temporary(ixc_);
  }

  /// End of solver
  void end(Block* block) throw();
  
  /// Access the Field message for buffering prolongation data
  FieldMsg ** pmsg_prolong(Block * block)
  {
    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();
    ScalarDescr *        scalar_descr = cello::scalar_descr_void();
    return (FieldMsg **)scalar_data->value(scalar_descr,i_msg_prolong_);
  }
  
  /// Access the Field message for buffering restriction data
  FieldMsg ** pmsg_restrict(Block * block, int ic)
  {
    ScalarData<void *> * scalar_data = block->data()->scalar_data_void();
    ScalarDescr *        scalar_descr = cello::scalar_descr_void();
    return (FieldMsg **)scalar_data->value(scalar_descr,i_msg_restrict_[ic]);
  }
  
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

protected: // methods

  FieldMsg * pack_field_
  (EnzoBlock *, int index_field, int refresh_type, int ic3[3]);

  void unpack_field_
  (EnzoBlock *, FieldMsg *, int index_field, int refresh_type);

  void copy_xc_to_x_(EnzoBlock *) throw();

protected: // attributes

  /// Matrix
  std::shared_ptr<Matrix> A_;

  /// Indices for coarse solver, domain solver, and smoother, prolong, restrict

  int index_solve_coarse_;
  int index_solve_domain_;
  int index_solve_smooth_;
  int index_prolong_;
  int index_restrict_;

  /// MG scalar id's
  int i_sync_restrict_;
  int i_sync_prolong_;
  int i_msg_restrict_[8];
  int i_msg_prolong_;

  /// Solver-specific temporary fields
  int ixc_;
  
  /// Block field attributes
  int mx_,my_,mz_;
  int gx_,gy_,gz_;

  /// The level of the coarse grid solve
  int coarse_level_;
};

#endif /* ENZO_ENZO_SOLVER_GRAVITY_DD_HPP */
