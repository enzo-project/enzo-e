// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_SimulationMpi.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Simulation] Declaration of the SimulationMpi class

#ifndef SIMULATION_SIMULATION_MPI_HPP
#define SIMULATION_SIMULATION_MPI_HPP

#if defined(CONFIG_USE_MPI) || ! defined(CONFIG_USE_CHARM)

class Block;
class Patch;

class SimulationMpi : public Simulation {

  /// @class    SimulationMpi
  /// @ingroup  Simulation
  /// @brief    [\ref Simulation] Simulation class for MPI

public: // functions

  /// Constructor
  SimulationMpi
  ( const char * parameter_file,
    const GroupProcess * group_process) throw();

  /// Destructor
  ~SimulationMpi() throw();

  /// Initialize the MPI Simulation
  virtual void initialize() throw();
  
  /// Run the simulation
  virtual void run() throw();

protected:

  void update_boundary_ (Block * block, bool boundary[3][2]) throw();
  void refresh_ghost_   (Block * block, Patch * patch, bool boundary[3][2]) throw();
  // void is_block_on_boundary_ (Block * block, bool boundary[3][2]) throw();
  
};

#endif /* ! CONFIG_USE_CHARM */

#endif /* SIMULATION_SIMULATION_MPI_HPP */

