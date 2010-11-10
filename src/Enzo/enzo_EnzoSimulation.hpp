// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_SIMULATION_HPP
#define ENZO_ENZO_SIMULATION_HPP

/// @file     enzo_EnzoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Declaration of the EnzoSimulation class

class EnzoSimulation : public Simulation {

  /// @class    EnzoSimulation
  /// @ingroup  Enzo
  /// @brief    Top-level description of user-implemented components

public: // interface

  /// Constructor
  EnzoSimulation(Global * global) throw();

  ~EnzoSimulation() throw()
  { delete enzo_; }

  /// Return the Enzo object created in EnzoSimulation's constructor
  EnzoDescr * enzo() throw ()
  { return enzo_; };

protected: // functions

  /// Read user parameters and initialize user objects
  void user_initialize_() throw ();

  /// Create named control method.
  MethodControl * create_method_control_ (std::string name_method_control) throw ();

  /// Create named timestep method.
  MethodTimestep * create_method_timestep_ (std::string name_method_timestep) throw ();

  /// Create named method method.
  MethodHyperbolic * create_method_hyperbolic_ (std::string name_method_hyperbolic) throw ();

private: // attributes

  EnzoDescr * enzo_;

};

#endif /* ENZO_ENZO_SIMULATION_HPP */

