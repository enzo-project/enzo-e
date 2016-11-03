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
  ( const char parameter_file[], int n);

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

public:

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

};

#endif /* ENZO_ENZO_SIMULATION_CHARM_HPP */
