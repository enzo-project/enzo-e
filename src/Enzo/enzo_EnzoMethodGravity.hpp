// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGravity.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2016-11-07
/// @brief    [\ref Enzo] Declaration of EnzoMethodGravity
///
/// Solve for gravitational potential using the specified linear solver

#ifndef ENZO_ENZO_METHOD_GRAVITY_HPP
#define ENZO_ENZO_METHOD_GRAVITY_HPP

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
  EnzoMethodGravity(const FieldDescr * field_descr,
		    int index_solver,
		    double grav_const);

  EnzoMethodGravity()
    : index_solver_(-1),
      grav_const_(0.0)
  {};

  /// Destructor
  ~EnzoMethodGravity() throw();

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGravity);
  
  /// Charm++ PUP::able migration constructor
  EnzoMethodGravity (CkMigrateMessage *m)
    : index_solver_(-1),
      grav_const_(0.0)
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

  }

  /// Solve for the gravitational potential
  virtual void compute( Block * block) throw();

  virtual std::string name () throw () 
  { return "gravity"; }

protected: // methods

  template <class T>
  void compute_ (EnzoBlock * enzo_block) throw();

protected: // attributes

  /// Solver index for the linear solver used to compute the potential
  int index_solver_;

  /// Gas constant, e.g. 6.67384e-8 (cgs)
  double grav_const_;
  
};


#endif /* ENZO_ENZO_METHOD_GRAVITY_HPP */
