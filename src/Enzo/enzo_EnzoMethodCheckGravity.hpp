// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_MethodCheckGravity.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     2020-01-21
/// @brief    [\ref Enzo] Declaration for the MethodCheckGravity class

#ifndef ENZO_ENZO_METHOD_CHECK_GRAVITY_HPP
#define ENZO_ENZO_METHOD_CHECK_GRAVITY_HPP

class EnzoMethodCheckGravity : public Method
{
  /// @class    EnzoMethodCheckGravity
  /// @ingroup  EnzoMethodCheckGravity
  /// @brief    [\ref EnzoMethodCheckGravity] Test method for basic Particle methods

public: // interface

  /// Create a new EnzoMethodCheckGravity
  EnzoMethodCheckGravity (std::string particle_type) throw() ;

  /// Destructor
  virtual ~EnzoMethodCheckGravity() throw()
  {};

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodCheckGravity);

  /// Charm++ PUP::able migration constructor
  EnzoMethodCheckGravity (CkMigrateMessage *m)
    : Method(m),
      particle_type_()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;
    Method::pup(p);
  };

public: // virtual functions

  /// Apply the method to advance a block one timestep 
  virtual void compute ( Block * block) throw();

  /// Return the name of this EnzoMethodCheckGravity
  virtual std::string name () throw ()
  { return "check_gravity"; }

protected: // functions


protected: // attributes

  std::string particle_type_;

};

#endif /* ENZO_ENZO_METHOD_CHECK_GRAVITY_HPP */
