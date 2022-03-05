// See LICENSE_CELLO file for license and copyright information


/// @file     enzo_EnzoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Enzo] Declaration of the EnzoSimulation class

#ifndef ENZO_ENZO_SIMULATION_CHARM_HPP
#define ENZO_ENZO_SIMULATION_CHARM_HPP

#include "charm++.h"
#include "enzo.decl.h"

class EnzoSimulation : public CBase_EnzoSimulation
			    
{

  /// @class    EnzoSimulation
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Simulation class for CHARM++ Enzo-E

public: // functions

  /// CHARM++ Constructor
  EnzoSimulation
  ( const char parameter_file[], int n);

  /// CHARM++ Constructor
  EnzoSimulation() : CBase_EnzoSimulation() {}

  /// CHARM++ Migration constructor
  EnzoSimulation(CkMigrateMessage * m) : CBase_EnzoSimulation(m)  
  {
  };

  /// Destructor
  virtual ~EnzoSimulation();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

#ifdef BUG_FIX_150
  /// Request by newly created EnzoBlock to get its MsgRefine object
  virtual void p_get_msg_refine(Index index);
#endif

  /// Barrier after constructor to ensure all EnzoSimulation objects created
  void r_startup_begun (CkReductionMsg *);

  /// EnzoMethodCheck
  void r_method_check_enter (CkReductionMsg *);
  void p_writer_created(int num_files, std::string ordering);
  
public: // virtual functions

  /// Initialize the Enzo Simulation
  virtual void initialize() throw();

  /// Return an EnzoFactory object, creating it if needed
  virtual const Factory * factory() const throw();

private: // functions

  virtual void initialize_config_() throw();

private: // attributes

  /// Checkpoint synchronization (depreciated)
  Sync sync_check_writer_created_;
  Sync sync_check_done_;
};

#endif /* ENZO_ENZO_SIMULATION_CHARM_HPP */
