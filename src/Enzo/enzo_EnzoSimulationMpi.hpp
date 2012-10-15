// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSimulationMpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Enzo] Declaration of the EnzoSimulationMpi class

#ifndef ENZO_ENZO_SIMULATION_MPI_HPP
#define ENZO_ENZO_SIMULATION_MPI_HPP

#if defined(CONFIG_USE_MPI) || ! defined(CONFIG_USE_CHARM)

class EnzoSimulationMpi : public SimulationMpi {

  /// @class    EnzoSimulationMpi
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for MPI Enzo-P

public: // functions

  /// Constructor
  EnzoSimulationMpi
  ( const char * parameter_file,
    const GroupProcess * group_process) throw();

  /// Destructor
  ~EnzoSimulationMpi() throw();

  /// Initialize the Enzo Simulation
  virtual void initialize() throw();

  /// Return an EnzoFactory object, creating it if needed
  virtual const Factory * factory() const throw();

private: // functions

  virtual void initialize_config_() throw();

};

#endif /* ! CONFIG_USE_CHARM */

#endif /* ENZO_ENZO_SIMULATION_MPI_HPP */

