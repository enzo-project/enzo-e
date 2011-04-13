// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_SIMULATION_MPI_HPP
#define ENZO_ENZO_SIMULATION_MPI_HPP

#ifndef CONFIG_USE_CHARM

/// @file     enzo_EnzoSimulationMpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Enzo] Declaration of the EnzoSimulationMpi class

class EnzoSimulationMpi : public EnzoSimulation {

  /// @class    EnzoSimulationMpi
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for MPI Enzo-P

public: // functions

  /// Constructor
  EnzoSimulationMpi
  ( const char * parameter_file,
    GroupProcess * group_process,
    int index) throw();

  /// Destructor
  ~EnzoSimulationMpi() throw();

  /// Run the simulation
  virtual void run() throw();

protected:
  
  GroupProcess * group_process_;

};

#endif /* ! CONFIG_USE_CHARM */

#endif /* ENZO_ENZO_SIMULATION_MPI_HPP */

