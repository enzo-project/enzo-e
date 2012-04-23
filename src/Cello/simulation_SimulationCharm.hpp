// See LICENSE_CELLO file for license and copyright information


/// @file     simulation_SimulationCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Simulation] Declaration of the SimulationCharm class

#ifndef SIMULATION_SIMULATION_CHARM_HPP
#define SIMULATION_SIMULATION_CHARM_HPP

#ifdef CONFIG_USE_CHARM

class SimulationCharm : public Simulation
			    
{

  /// @class    SimulationCharm
  /// @ingroup  Simulation
  /// @brief    [\ref Simulation] Simulation class for CHARM++

public: // functions

  /// CHARM++ Constructor
  SimulationCharm
  ( const char parameter_file[], int n) throw();

  /// Destructor
  ~SimulationCharm() throw();

  /// Initialize the Enzo Simulation
  virtual void initialize() throw();

  /// Run the simulation
  virtual void run() throw();

};

#endif /* CONFIG_USE_CHARM */

#endif /* SIMULATION_SIMULATION_CHARM_HPP */
