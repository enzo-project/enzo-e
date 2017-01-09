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
   std::string smooth,
   double      smooth_weight,
   int         smooth_pre,
   int         smooth_coarse,
   int         smooth_post,
   bool is_singular,
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
       smooth_(""),
       smooth_pre_(NULL),
       smooth_coarse_(NULL),
       smooth_post_(NULL),
       restrict_(NULL),
       prolong_(NULL),
       is_singular_(false),
       rank_(0),
       iter_max_(0), 
       monitor_iter_(0),
       ib_(0), ic_(0), id_(0), ir_(0), ix_(0),
       min_level_(0),
       max_level_(0),
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
    p | smooth_;
    p | smooth_pre_;
    p | smooth_coarse_;
    p | smooth_post_;
    p | restrict_;
    p | prolong_;
    p | is_singular_;
    p | rank_;
    p | iter_max_;
    p | monitor_iter_;

    p | ib_;
    p | ic_;
    p | id_;
    p | ir_;
    p | ix_;

    p | min_level_;
    p | max_level_;

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
  void restrict_recv(EnzoBlock * enzo_block) throw();

  /// Prolong the correction C to the next-finer level
  template <class T>
  void prolong_recv(EnzoBlock * enzo_block) throw();

  /// Apply post-smoothing to the current level
  template <class T>
  void post_smooth(EnzoBlock * enzo_block) throw();

  void set_bs(long double bs) throw() { bs_ = bs; }
  void set_bc(long double bc) throw() { bc_ = bc; }

  template <class T>
  void begin_solve_(EnzoBlock * enzo_block) throw();

protected: // methods

  template <class T>
  void enter_solver_(EnzoBlock * enzo_block) throw();

  template <class T>
  void begin_cycle_(EnzoBlock * enzo_block) throw();
  
  /// Solve the coarse-grid equation A*C = R
  template <class T>
  void solve_coarse_(EnzoBlock * enzo_block) throw();

  /// Prolong the correction C to the next-finer level
  template <class T>
  void prolong_send_(EnzoBlock * enzo_block) throw();

  template <class T>
  void end_cycle_(EnzoBlock * enzo_block) throw();
  
  void monitor_output_(EnzoBlock * enzo_block) throw();

  bool is_converged_(EnzoBlock * enzo_block) const;

protected: // attributes

  /// Matrix
  Matrix * A_;

  /// Smoother name
  std::string smooth_;
  
  /// Pre-smoother
  Solver * smooth_pre_;

  /// Coarse smoother
  Solver * smooth_coarse_;

  /// Post smoother
  Solver * smooth_post_;

  /// Restriction
  Restrict * restrict_;

  /// Prolongation
  Prolong * prolong_;

  /// Whether you need to subtract of the nullspace of A from b, e.g. fully
  /// periodic or Neumann problems
  bool is_singular_;

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
  int id_;
  int ir_;
  int ix_;

  /// Minimum refinement level (may be < 0)
  int min_level_;

  /// Maximum refinement level
  int max_level_;

  /// Block field attributes
  int mx_,my_,mz_;
  int nx_,ny_,nz_;
  int gx_,gy_,gz_;

  /// scalars used for projections of singular systems
  long double bs_;
  long double bc_;
};

#endif /* ENZO_ENZO_SOLVER_GRAVITY_MG0_HPP */
