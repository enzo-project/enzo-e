// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravityMg0.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2015-06-02
/// @brief    [\ref Enzo] Declaration of EnzoMethodGravityMg0
///
/// Multigrid for solving a linear system on the root-level grid only

#ifndef ENZO_ENZO_METHOD_GRAVITY_MG0_HPP
#define ENZO_ENZO_METHOD_GRAVITY_MG0_HPP

class EnzoMethodGravityMg0 : public Method {

  /// @class    EnzoMethodGravityMg0
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Multigrid on the root-level grid.  For use either
  /// as a Gravity solver on non-adaptive problems, or as a preconditioner
  /// for a Krylov subspace method, as in Dan Reynold's HG method.

public: // interface

  /// Create a new EnzoMethodGravityMg0 object
  EnzoMethodGravityMg0
  (const FieldDescr * field_descr, int rank,
   double grav_const,
   int iter_max, 
   int monitor_iter,
   bool is_singular,
   Compute * smooth,
   Restrict * restrict,
   Prolong * prolong,
   int min_level,
   int max_level);

  EnzoMethodGravityMg0() {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravityMg0);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravityMg0 (CkMigrateMessage *m) {}

  /// Destructor
  ~EnzoMethodGravityMg0 () throw();

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
    p | monitor_iter_;
    p | rr_;
    p | rr0_;

    p | irho_;
    p | iphi_;
    p | ib_;
    p | ir_;
    p | ix_;
    p | ic_;

    p | min_level_;
    p | max_level_;

    p | mx_;
    p | my_;
    p | mz_;
  }

  /// Solve for the gravitational potential
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "gravity_mg"; }

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
  template <class T>
  void exit_solver_(EnzoBlock * enzo_block, int retval) throw();
  
  void monitor_output_(EnzoBlock * enzo_block) throw();

  bool is_converged_(EnzoBlock * enzo_block) const;

  template <class T>
  void zaxpy_ (T * Z, double a, const T * X, const T * Y) const throw();


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

  /// Maximum number of MG iterations
  int iter_max_;

  /// How often to display progress
  int monitor_iter_;

  /// Current residual
  long double rr_;

  /// Initial residual
  long double rr0_;

  /// Density and potential field id's

  int irho_;
  int iphi_;

  /// MG vector id's
  int ib_;
  int ir_;
  int ix_;
  int ic_;

  /// Minimum refinement level (may be < 0)
  int min_level_;

  /// Maximum refinement level
  int max_level_;

  /// Block field attributes
  int mx_,my_,mz_;
};

#endif /* ENZO_ENZO_METHOD_GRAVITY_MG0_HPP */
