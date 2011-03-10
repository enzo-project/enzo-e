// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_SIMULATION_SERIAL_HPP
#define ENZO_ENZO_SIMULATION_SERIAL_HPP

/// @file     enzo_EnzoSimulationSerial.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Enzo] Declaration of the EnzoSimulationSerial class

class EnzoSimulationSerial : public Simulation {

  /// @class    EnzoSimulationSerial
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for Enzo

public: // interface

  /// Constructor
  EnzoSimulationSerial(Parameters * parameters) throw();

  /// Destructor
  ~EnzoSimulationSerial() throw();

  /// Override Simulation initialize
  virtual void initialize() throw ();

  /// Finalize the Simulation after running it
  virtual void finalize() throw();

  /// Run the simulation
  virtual void run() throw();

  /// Load a Simulation from disk
  virtual void read() throw();

  /// Write a Simulation state to disk
  virtual void write() throw();

  /// Create a new Mesh: FACTORY METHOD DESIGN PATTERN
  virtual Mesh * create_mesh (int nx,int ny,int nz,
			      int nbx,int nby,int nbz) throw()
  { 
    return new EnzoMesh (nx,ny,nz,nbx,nby,nbz);
  };

public: // functions

  

protected: // virtual functions

  /// Create named stopping object
  Stopping * create_stopping_ (std::string name) throw ();

  /// Create named timestep object
  Timestep * create_timestep_ (std::string name) throw ();

  /// Create named initial conditions object
  Initial * create_initial_ (std::string name) throw ();

  /// Create named boundary conditions object
  Boundary * create_boundary_ (std::string name) throw ();

  /// Create named method object
  Method * create_method_ (std::string name) throw ();

private: // functions

  /// Output data
  void output_images_( EnzoBlock * enzo_block,
		       const char * file_format,
		       int cycle,
		       int cycle_skip=1) throw ();

  void deallocate_() throw();

};

#endif /* ENZO_ENZO_SIMULATION_SERIAL_HPP */

