// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_SIMULATION_HPP
#define ENZO_ENZO_SIMULATION_HPP

/// @file     enzo_EnzoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Enzo] Declaration of the EnzoSimulation class

class EnzoSimulation : public Simulation {

  /// @class    EnzoSimulation
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for Enzo

public: // interface

  /// Constructor
  EnzoSimulation(Error   * error,
		 Monitor * monitor) throw();

  /// Destructor
  ~EnzoSimulation() throw()
  { delete enzo_descr_; }

  /// Return the Enzo object created in EnzoSimulation's constructor
  EnzoDescr * enzo() throw ()
  { return enzo_descr_; };

  /// Override Simulation initialize
  void initialize(FILE * fp) throw ();

protected: // virtual functions

  /// Create named control object
  Control * create_control_ (std::string name) throw ();

  /// Create named timestep object
  Timestep * create_timestep_ (std::string name) throw ();

  /// Create named initial conditions object
  Initial * create_initial_ (std::string name) throw ();

  /// Create named boundary conditions object
  Boundary * create_boundary_ (std::string name) throw ();

  /// Create named method object
  Method * create_method_ (std::string name) throw ();

private: // attributes

  EnzoDescr * enzo_descr_;

};

#endif /* ENZO_ENZO_SIMULATION_HPP */

