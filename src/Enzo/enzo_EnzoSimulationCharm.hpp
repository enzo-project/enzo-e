// $Id: enzo_EnzoSimulationCharm.hpp 2115 2011-03-17 19:59:55Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_SIMULATION_CHARM_HPP
#define ENZO_ENZO_SIMULATION_CHARM_HPP

/// @file     enzo_EnzoSimulationCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Enzo] Declaration of the EnzoSimulationCharm class

class EnzoSimulationCharm : public EnzoSimulation {

  /// @class    EnzoSimulationCharm
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for CHARM++ Enzo-P

public: // functions

  /// Constructor
  EnzoSimulationCharm
  (
   Parameters * parameters,
   GroupProcess *
   ) throw();

  /// Destructor
  ~EnzoSimulationCharm() throw();

  /// Run the simulation
  virtual void run() throw();
};

#endif /* ENZO_ENZO_SIMULATION_CHARM_HPP */

