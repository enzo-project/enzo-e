// See LICENSE_CELLO file for license and copyright information


/// @file     simulation_SimulationCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Simulation] Declaration of the SimulationCharm class

#ifndef SIMULATION_SIMULATION_CHARM_HPP
#define SIMULATION_SIMULATION_CHARM_HPP

class SimulationCharm : public Simulation
			    
{

  /// @class    SimulationCharm
  /// @ingroup  Simulation
  /// @brief    [\ref Simulation] Simulation class for CHARM++

public: // functions

  /// CHARM++ Constructor
  SimulationCharm
  ( const char parameter_file[], int n) throw();

  /// Destructor
  ~SimulationCharm() throw();

  /// CHARM++ Constructor
  SimulationCharm() 
  {
  }

  /// CHARM++ Migration constructor
  SimulationCharm(CkMigrateMessage*m)
    : Simulation(m)
  {}

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    Simulation::pup(p);
    p | sync_output_begin_;
    p | sync_output_write_;

  // clear sync if unpacking: load balancing expects syncs to be
  // updated by CommBlock(CkMigrateMessage) (increment) and
  // ~CommBlock() (decrement), but checkpoint / restart then
  // double-counts Blocks.

    const bool up = p.isUnpacking();

    if (up) sync_output_begin_.set_stop(0);
    if (up) sync_output_write_.set_stop(0);
  }

  /// Initialize the Charm++ Simulation
  virtual void initialize() throw();

  // /// Run the simulation
  // virtual void run() throw();

  /// Add a new CommBlock to this local branch
  void insert_block() ;

  /// Remove a CommBlock from this local branch
  void delete_block() ;

  /// Wait for all Hierarchy to be initialized before creating any CommBlocks
  void r_initialize_forest(CkReductionMsg * msg);

  /// Wait for all local patches to be created before calling run
  void r_initialize_hierarchy(CkReductionMsg * msg);

  /// Call output on Problem list of Output objects
  void p_begin_output()
  { begin_output(); }
  void begin_output ();
  void output_exit();
  void r_output(CkReductionMsg * msg);

  //  void r_output (CkReductionMsg * msg);

  /// Reduce output, using p_output_write to send data to writing processes
  void s_write() { write_(); };
  void write_();

  /// Continue on to Problem::output_wait()
  void r_write(CkReductionMsg * msg);

  /// Continue on to Problem::output_wait() from checkpoint
  void r_write_checkpoint();

  /// Receive data from non-writing process, write to disk, close, and
  /// proceed with next output
  void p_output_write (int n, char * buffer);

  void compute ();

  void p_monitor();

  void p_monitor_performance()
  { monitor_performance(); };

  virtual void monitor_performance();

  /// Reduction for performance data
  void r_monitor_performance (CkReductionMsg * msg);

 
protected: // attributes

  Sync sync_output_begin_;
  Sync sync_output_write_;

#ifdef CONFIG_USE_MEMORY
  int trace_mem_;
#endif

};

#endif /* SIMULATION_SIMULATION_CHARM_HPP */
