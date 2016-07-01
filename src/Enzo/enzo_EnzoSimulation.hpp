// See LICENSE_CELLO file for license and copyright information


/// @file     enzo_EnzoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Enzo] Declaration of the EnzoSimulation class

#ifndef ENZO_ENZO_SIMULATION_CHARM_HPP
#define ENZO_ENZO_SIMULATION_CHARM_HPP

#include "charm++.h"
#include "enzo.decl.h"

class EnzoSimulation : public BASE_ENZO_SIMULATION
			    
{

  /// @class    EnzoSimulation
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for CHARM++ Enzo-P

public: // functions

  /// CHARM++ Constructor
  EnzoSimulation
  ( const char parameter_file[], int n, int node_size);

  /// CHARM++ Constructor
  EnzoSimulation() : BASE_ENZO_SIMULATION() {}

  /// CHARM++ Migration constructor
  EnzoSimulation(CkMigrateMessage * m) : BASE_ENZO_SIMULATION(m)  
  {
  };

  /// Destructor
  ~EnzoSimulation();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Barrier after constructor to ensure all EnzoSimulation objects created
  void r_startup_begun (CkReductionMsg *);

  /// Read parameter file one process at a time
  void p_read_parameters()
  { 
    read_parameters_(); 
  }
private:
  void read_parameters_();
public:

  /// Barrier after parameters read to synchronize before initializing 
  void r_startup_finished (CkReductionMsg *);

  /// Continue on to Problem::output_wait() from checkpoint
  virtual void r_write_checkpoint();

public: // virtual functions

  /// Initialize the Enzo Simulation
  virtual void initialize() throw();

  /// Return an EnzoFactory object, creating it if needed
  virtual const Factory * factory() const throw();

private: // functions

  virtual void initialize_config_() throw();

private: // attributes

  // Number of processes per node: default is MAX_NODE_SIZE, but can be
  // overridden using NODE_SIZE file.  Used to serialize reading parameters
  // within a node.
  int node_size_; 

};

#endif /* ENZO_ENZO_SIMULATION_CHARM_HPP */
