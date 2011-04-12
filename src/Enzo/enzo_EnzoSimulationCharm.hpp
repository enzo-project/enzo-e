// $Id: enzo_EnzoSimulationCharm.hpp 2115 2011-03-17 19:59:55Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_SIMULATION_CHARM_HPP
#define ENZO_ENZO_SIMULATION_CHARM_HPP

/// @file     enzo_EnzoSimulationCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Enzo] Declaration of the EnzoSimulationCharm class

#ifdef CONFIG_USE_CHARM

#include PARALLEL_CHARM_INCLUDE(enzo.decl.h)

class EnzoSimulationCharm : public CBase_EnzoSimulationCharm,
			    public EnzoSimulation
{

  /// @class    EnzoSimulationCharm
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for CHARM++ Enzo-P

public: // functions

  /// CHARM++ Constructor
  EnzoSimulationCharm
  ( const char parameter_file[], int n, int index_simulation) throw();

  /// Constructor
  EnzoSimulationCharm
  ( const char * parameter_file,
    GroupProcess * ) throw();

  /// Destructor
  ~EnzoSimulationCharm() throw();

  /// Run the simulation
  virtual void run() throw();

  //==================================================

  /// Apply Output and Monitor before each cycle
  void p_prepare() throw();

  /// Request all Mesh blocks to send output to main::p_output_close()
  void p_output(int index, int cycle, double time) throw();

private:
  int count_prepare_;

};

#endif /* CONFIG_USE_CHARM */

#endif /* ENZO_ENZO_SIMULATION_CHARM_HPP */
