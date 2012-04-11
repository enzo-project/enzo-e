// See LICENSE_CELLO file for license and copyright information


/// @file     enzo_EnzoSimulationCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Enzo] Declaration of the EnzoSimulationCharm class

#ifndef ENZO_ENZO_SIMULATION_CHARM_HPP
#define ENZO_ENZO_SIMULATION_CHARM_HPP

#ifdef CONFIG_USE_CHARM

class EnzoSimulationCharm : public EnzoSimulation
			    
{

  /// @class    EnzoSimulationCharm
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for CHARM++ Enzo-P

public: // functions

  /// CHARM++ Constructor
  EnzoSimulationCharm
  ( const char parameter_file[], int n) throw();

  /// Destructor
  ~EnzoSimulationCharm() throw();

  /// Run the simulation
  virtual void run() throw();

};

#endif /* CONFIG_USE_CHARM */

#endif /* ENZO_ENZO_SIMULATION_CHARM_HPP */
