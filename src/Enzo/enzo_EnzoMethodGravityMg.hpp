// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityMg.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2015-04-24 22:32:18
/// @brief    [\ref Enzo] Declaration of EnzoMethodGravityMg
///
/// Multigrid on adaptive meshes for solving for self-gravity on field
/// data.

#ifndef ENZO_ENZO_METHOD_GRAVITY_MG_HPP
#define ENZO_ENZO_METHOD_GRAVITY_MG_HPP

class EnzoMethodGravityMg : public Method {

  /// @class    EnzoMethodGravityMg
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Demonstration method to solve self-gravity
  /// using Multigrid.

public: // interface

  /// Create a new EnzoMethodGravityMg object
  EnzoMethodGravityMg(FieldDescr * field_descr, int rank,
		      double grav_const,
		      int iter_max, 
		      double res_tol,
		      int monitor_iter,
		      bool is_singular,
		      Compute * smooth,
		      Restrict * restrict,
		      Prolong * prolong,
		      int level_min,
		      int level_max);

  EnzoMethodGravityMg() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityMg);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityMg (CkMigrateMessage *m) {}

  /// Destructor
  ~EnzoMethodGravityMg () throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    // NOTE: change this function whenever attributes change

    TRACEPUP;

    Method::pup(p);

    p | A_;
    p | smooth_;
    p | restrict_;
    p | prolong_;
    p | is_singular_;
    p | rank_;
    p | grav_const_;
    p | iter_max_;
    p | res_tol_;
    p | monitor_iter_;
    p | rr_;
    p | rr0_;
    p | rr_min_;
    p | rr_max_;

    p | idensity_;
    p | ipotential_;
    p | ib_;
    p | ir_;
    p | ix_;
    p | iy_;

    p | nx_;
    p | ny_;
    p | nz_;

    p | mx_;
    p | my_;
    p | mz_;

    p | gx_;
    p | gy_;
    p | gz_;

    p | iter_;
    p | is_active_;
    p | level_;
    p | level_min_;
    p | level_max_;
  }

  /// Solve for the gravitational potential
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "gravity_mg"; }

protected: // methods

  /// Apply pre-smoothing on the current level
  void pre_smooth (EnzoBlock * enzo_block, int level);
  /// Compute the residual R = B - A*X on the current level
  void compute_residual (EnzoBlock * enzo_block, int level);
  /// Restrict the residual R to the next-coarser level
  void restrict_residual (EnzoBlock * enzo_block, int level);
  /// Solve the coarse-grid equation A*C = R
  void coarse_solve (EnzoBlock * enzo_block, int level);
  /// Prolong the correction C to the next-finer level
  void prolong_correction (EnzoBlock * enzo_block, int level);
  /// Update the solution X <= X + C on the current level
  void update_solution (EnzoBlock * enzo_block, int level);
  /// Apply post-smoothing to the current level
  void post_smooth (EnzoBlock * enzo_block, int level);

  void monitor_output_(EnzoBlock * enzo_block) throw();

  template <class T>
  void compute_ (EnzoBlock * enzo_block) throw();

  template <class T>
  void mg_end (EnzoBlock * enzo_block, int retval) throw();

  /// Set whether current Block is a leaf--if not don't touch data
  void set_active(Block * block) throw()
  { is_active_ = (level_ == block->level()); }

  bool is_active () const throw()
  { return is_active_; }

protected: // attributes

  /// Matrix
  Matrix * A_;

  /// Smoother
  Compute * smooth_;

  /// Restriction
  Restrict * restrict_;

  /// Prolongation
  Prolong * prolong_;

  /// Whether you need to subtract of the nullspace of A from b, e.g. fully
  /// periodic or Neumann problems
  bool is_singular_;

  /// Dimensionality of the problem
  int rank_;

  /// Gas constant, e.g. 6.67384e-8 (cgs)
  double grav_const_;

  /// Maximum number of Mg iterations
  int iter_max_;

  /// Convergence tolerance on the residual reduction rz_ / rz0_
  double res_tol_;

  /// How often to display progress
  int monitor_iter_;

  /// Current residual
  long double rr_;

  /// Initial residual
  long double rr0_;

  /// Minimum residual
  long double rr_min_;

  /// Maximum residual
  long double rr_max_;

  /// Density and potential field id's

  int idensity_;
  int ipotential_;

  /// MG vector id's
  int ib_;
  int ir_;
  int ix_;
  int iy_;

  /// Block field attributes
  int nx_,ny_,nz_;
  int mx_,my_,mz_;
  int gx_,gy_,gz_;

  /// Current MG iteration
  int iter_;

  /// Whether current block is active (set using set_active)
  int is_active_;

  /// Current active level
  int level_;

  /// Minimum refinement level (may be < 0)
  int level_min_;

  /// Maximum refinement level
  int level_max_;

};

#endif /* ENZO_ENZO_METHOD_GRAVITY_MG_HPP */
