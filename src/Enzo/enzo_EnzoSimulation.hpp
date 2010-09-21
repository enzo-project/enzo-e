// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_USER_DESCR_HPP
#define ENZO_ENZO_USER_DESCR_HPP

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
  UserControl * create_user_control_ (std::string name_user_control) throw ();

  /// Create named timestep method.
  UserTimestep * create_user_timestep_ (std::string name_user_timestep) throw ();

  /// Create named user method.
  UserMethod * create_user_method_ (std::string name_user_method) throw ();

private: // attributes

  EnzoDescr * enzo_;

};

#endif /* ENZO_ENZO_USER_DESCR_HPP */

