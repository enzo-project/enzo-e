// See LICENSE_CELLO file for license and copyright information

#ifndef IO_IO_SIMULATION_HPP
#define IO_IO_SIMULATION_HPP

/// @file     io_IoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Aug 29 16:53:50 PDT 2011
/// @brief    [\ref Io] Implement disk IO for simulations

class Simulation;

class IoSimulation {

  /// @class    IoSimulation
  /// @ingroup  Io
  /// @brief [\ref Io] IoSimulation is used for writing and
  /// reading Simulation data to and from disk

public: // interface

  /// Constructor
  IoSimulation(Simulation * simulation) throw();

  /// Read a simulation from disk
  void read () throw();

  /// Write a Simulation to disk
  void write () throw();

private: // functions

private: // attributes

  /// Pointer to the simulation
  Simulation * simulation_; 


};

#endif /* IO_IO_SIMULATION_HPP */

