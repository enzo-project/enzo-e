// See LICENSE_CELLO file for license and copyright information


/// @file     enzo_EnzoSimulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Enzo] Declaration of the EnzoSimulation class

#ifndef ENZO_ENZO_SIMULATION_CHARM_HPP
#define ENZO_ENZO_SIMULATION_CHARM_HPP

class CProxy_IoEnzoReader;
class CProxy_IoEnzoWriter;

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

#ifdef BYPASS_CHARM_MEM_LEAK
  /// Request by newly created EnzoBlock to get its MsgRefine object
  virtual void p_get_msg_refine(Index index);

  virtual void p_get_msg_check(Index index);
  void set_msg_check (Index index, EnzoMsgCheck *);
  EnzoMsgCheck * get_msg_check (Index index);
#endif

  /// Barrier after constructor to ensure all EnzoSimulation objects created
  void r_startup_begun (CkReductionMsg *);

  /// EnzoMethodCheck
  void r_method_check_enter (CkReductionMsg *);
  void p_check_done();
  void p_set_io_reader(CProxy_IoEnzoReader io_reader);
  void p_set_io_writer(CProxy_IoEnzoWriter io_writer);
  void set_sync_check_writer(int count)
  { sync_check_writer_created_.set_stop(count); }
  void p_io_reader_created();

  /// Read in and initialize the next refinement level from a checkpoint;
  /// or exit if done
  void p_restart_next_level();
  void p_restart_level_created();

public: // virtual functions

  /// Initialize the Enzo Simulation
  virtual void initialize() throw();

  /// Return an EnzoFactory object, creating it if needed
  virtual const Factory * factory() const throw();

private: // functions


private: // virtual functions

  virtual void initialize_config_() throw();

private: // attributes

  /// Checkpoint synchronization
  Sync                     sync_check_writer_created_;
  Sync                     sync_check_done_;
  int                      check_num_files_;
  std::string              check_ordering_;
  std::vector<std::string> check_directory_;

  /// Current restart level
  int restart_level_; 
#ifdef BYPASS_CHARM_MEM_LEAK
  std::map<Index,EnzoMsgCheck *> msg_check_map_;
#endif
};

#endif /* ENZO_ENZO_SIMULATION_CHARM_HPP */
