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
  { delete enzo_; }

  /// Override Simulation initialize
  virtual void initialize(FILE * parameter_file) throw ();

  /// Finalize the Simulation after running it
  virtual void finalize() throw();

  /// Run the simulation
  virtual void run() throw();

  /// Load a Simulation from disk
  virtual void read() throw();

  /// Write a Simulation state to disk
  virtual void write() throw();


public: // functions

  /// Return the Enzo object created in EnzoSimulation's constructor
  EnzoDescr * enzo() throw ()
  { return enzo_; };

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

  EnzoDescr * enzo_;

};

#endif /* ENZO_ENZO_SIMULATION_HPP */

