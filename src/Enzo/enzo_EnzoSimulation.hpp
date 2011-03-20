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
  /// @brief    [\ref Enzo] Simulation base class for Enzo-P

public: // functions

  /// Constructor
  EnzoSimulation
  (
   Parameters * parameters,
   GroupProcess *
   ) throw();

  /// Destructor
  virtual ~EnzoSimulation() throw();

  /// Override Simulation initialize
  virtual void initialize() throw ();

  /// Finalize the Simulation after running it
  virtual void finalize() throw();

  /// Run the simulation
  virtual void run() throw() = 0;

  /// Load a Simulation from disk
  virtual void read() throw();

  /// Write a Simulation state to disk
  virtual void write() const throw();

protected: // virtual functions

  /// Create named stopping object
  virtual Stopping * create_stopping_ (std::string name) throw ();

  /// Create named timestep object
  virtual Timestep * create_timestep_ (std::string name) throw ();

  /// Create named initial conditions object
  virtual Initial * create_initial_ (std::string name) throw ();

  /// Create named boundary conditions object
  virtual Boundary * create_boundary_ (std::string name) throw ();

  /// Create named method object
  virtual Method * create_method_ (std::string name) throw ();

  /// Create output object for the given filename
  virtual Output * create_output_ (std::string name) throw ();

};

#endif /* ENZO_ENZO_SIMULATION_HPP */

