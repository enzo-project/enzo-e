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

  /// Destructor
  ~EnzoSimulation() throw()
  { delete enzo_descr_; }

  /// Return the Enzo object created in EnzoSimulation's constructor
  EnzoDescr * enzo() throw ()
  { return enzo_descr_; };

  /// Override Simulation initialize
  void initialize(std::string parameter_file) throw ();

protected: // virtual functions

  /// Create named control method.
  MethodControl * create_control_ (std::string name_control) throw ();

  /// Create named timestep method.
  MethodTimestep * create_timestep_ (std::string name_timestep) throw ();

  /// Create named method method.
  MethodHyperbolic * create_method_ (std::string name_method) throw ();

  /// Create named method method.
  MethodInitial * create_initial_ (std::string name_initial) throw ();

private: // attributes

  EnzoDescr * enzo_descr_;

};

#endif /* ENZO_ENZO_SIMULATION_HPP */

