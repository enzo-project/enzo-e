// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2016-11-07
/// @brief    [\ref Enzo] Declaration of EnzoMethodGravity
///
/// Solve for gravitational potential using the specified linear solver

#ifndef ENZO_ENZO_METHOD_GRAVITY_HPP
#define ENZO_ENZO_METHOD_GRAVITY_HPP

enum TypeSuper {
  super_type_potential = 0,
  super_type_accelerations = 1 };

class EnzoMethodGravity : public Method {

  /// @class    EnzoMethodGravity
  /// @ingroup  Enzo
  ///
  /// @brief [\ref Enzo] Method to solve self-gravity.  Uses the
  /// "density_total" field, which must be initialized with the
  /// density field(s) and particles with "mass" attribute or
  /// constant.  Applies the solver to solve for the "potential"
  /// field.

public: // interface

  /// Create a new EnzoMethodGravity object
  EnzoMethodGravity(int index_solver,
		    double grav_const,
		    int order,
		    bool accumulate,
		    int index_prolong,
		    double dt_max,
                    std::string type_super);

  EnzoMethodGravity()
    : index_solver_(-1),
      grav_const_(0.0),
      order_(4),
      ir_exit_(-1),
      index_prolong_(0),
      dt_max_(0.0),
      ipotential_curr_(-1),
      ipotential_prev_(-1),
      is_time_curr_(-1),
      is_time_prev_(-1),
      type_super_(super_type_potential)
  {};

  /// Destructor
  virtual ~EnzoMethodGravity() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravity);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravity (CkMigrateMessage *m)
    : Method (m),
      index_solver_(-1),
      grav_const_(0.0),
      order_(4),
      ir_exit_(-1),
      index_prolong_(0),
      dt_max_(0.0),
      ipotential_curr_(-1),
      ipotential_prev_(-1),
      is_time_curr_(-1),
      is_time_prev_(-1),
      type_super_(super_type_potential)
  { }

  /// CHARM++ Pack / Unpack function
//----------------------------------------------------------------------

  void pup (PUP::er &p)
  {

    // NOTE: change this function whenever attributes change

    TRACEPUP;

    Method::pup(p);

    p | index_solver_;
    p | grav_const_;
    p | order_;
    p | dt_max_;
    p | ir_exit_;
    p | ipotential_curr_;
    p | ipotential_prev_;
    p | is_time_curr_;
    p | is_time_prev_;
    p | type_super_;

  }

  /// Solve for the gravitational potential
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "gravity"; }

  /// Compute maximum timestep for this method
  virtual double timestep (Block * block) throw() ;

  /// Compute accelerations from potential and exit solver
  void compute_accelerations (EnzoBlock * enzo_block) throw();

  void refresh_potential (EnzoBlock * enzo_block) throw();

  protected: // methods

  void compute_ (EnzoBlock * enzo_block) throw();

  /// Compute maximum timestep for this method
  double timestep_ (Block * block) throw() ;

  /// Scalar times for previous and current saved potentials when max_supercycling > 1
  double & s_time_curr_(Block * block)
  { return *block->data()->scalar_double().value(is_time_curr_); }
  double & s_time_prev_(Block * block)
  { return *block->data()->scalar_double().value(is_time_prev_); }

  void refresh_add_potentials_(Refresh *);
  void refresh_add_accelerations_(Refresh *);
  void super_save_potential_(EnzoBlock * );
  void super_save_accelerations_(EnzoBlock * );
  void super_extrapolate_potential_(EnzoBlock * );
  void super_extrapolate_accelerations_(EnzoBlock * );

protected: // attributes

  /// Solver index for the linear solver used to compute the potential
  int index_solver_;

  /// Gravity constant, e.g. 6.67384e-8 (cgs)
  double grav_const_;

  /// Order of Laplacian and acceleration computation: 2 or 4
  /// (Note EnzoMatrixLaplacian supports order=6 as well)
  int order_;

  /// Refresh id's
  int ir_exit_;

  /// Prolongation
  int index_prolong_;

  /// Maximum timestep
  double dt_max_;

  /// Temporary fields
  int ipotential_curr_;
  int ipotential_prev_;

  /// Scalars
  int is_time_curr_;
  int is_time_prev_;

  /// Type of super-cycling (if any): 0 potential 1:accelerations
  int type_super_;
};


#endif /* ENZO_ENZO_METHOD_GRAVITY_HPP */
