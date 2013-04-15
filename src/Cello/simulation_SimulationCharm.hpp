// See LICENSE_CELLO file for license and copyright information


/// @file     simulation_SimulationCharm.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-03-17
/// @brief    [\ref Simulation] Declaration of the SimulationCharm class

#ifndef SIMULATION_SIMULATION_CHARM_HPP
#define SIMULATION_SIMULATION_CHARM_HPP

#ifdef CONFIG_USE_CHARM

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
  SimulationCharm() {}

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
    p | block_loop_;
  }

  /// Initialize the Charm++ Simulation
  virtual void initialize() throw();

  /// Run the simulation
  virtual void run() throw();

  /// Add a new CommBlock to this local branch
  inline void insert_block() 
  {
    WARNING("SimulationCharm::insert_block()",
	    "Migrating CommBlocks will disturb local counts");
    ++block_loop_;
    TRACE2 ("++local count %d = %d",group_process_->rank(),block_loop_.stop());
  }

  /// Remove a CommBlock from this local branch
  inline void delete_block() 
  {
    WARNING("SimulationCharm::delete_block()",
	    "Migrating CommBlocks will disturb local counts");
    --block_loop_; 
    TRACE2 ("--local count %d = %d",group_process_->rank(),block_loop_.stop());
  }

  /// Call initialize()
  void p_initialize() { initialize(); }

  /// Wait for all local patches to be created before calling run
  void s_initialize();

  // Call initialization on Problem list of Initial objects
  void p_initial ();
  void initial ();

  // Call output on Problem list of Output objects
  void p_output ();

  /// Reduce output, using p_output_write to send data to writing processes
  void s_write();

  /// Receive data from non-writing process, write to disk, close, and
  /// proceed with next output
  void p_output_write (int n, char * buffer);

  // /// Wait for all local patches to check in before proceeding
  // void s_patch(CkCallback function);

  /// Wait for all local patches to check in before proceeding to refresh
  void s_initial();

  /// Refresh ghost zones (Charm++ entry)
  void p_refresh();

  /// Refresh ghost zones
  void refresh();

  // Stopping criteria and computation
  void c_compute ();

protected: // attributes
  Loop block_loop_;

};

#endif /* CONFIG_USE_CHARM */

#endif /* SIMULATION_SIMULATION_CHARM_HPP */
