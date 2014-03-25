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
    p | block_sync_;
  }

  /// Initialize the Charm++ Simulation
  virtual void initialize() throw();

  // /// Run the simulation
  // virtual void run() throw();

  /// Add a new CommBlock to this local branch
  inline void insert_block() 
  {
#ifdef CELLO_DEBUG
    PARALLEL_PRINTF ("%d: ++block_sync_ %d\n",  CkMyPe(),block_sync_.stop());
#endif
    ++block_sync_;
  }

  /// Remove a CommBlock from this local branch
  inline void delete_block() 
  {
#ifdef CELLO_DEBUG
    PARALLEL_PRINTF ("%d: --block_sync_ %d\n",  CkMyPe(),block_sync_.stop());
#endif
    --block_sync_; 
  }

  inline int block_count() const
  { 
#ifdef CELLO_DEBUG
    PARALLEL_PRINTF ("%d: block_sync_ %d\n",  CkMyPe(),block_sync_.stop());
#endif
    return block_sync_.stop(); 
  }

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

  /// Receive data from non-writing process, write to disk, close, and
  /// proceed with next output
  void p_output_write (int n, char * buffer);

  void compute ();


  /// Output Performance information to stdout
  void p_performance_output()
  { performance_output(); }

  virtual void performance_output();

  /// Reduction for performance data
  void r_performance_reduce (CkReductionMsg * msg);

 
protected: // attributes

  Sync block_sync_;

#ifdef TRACE_MEMORY
  int trace_mem_;
#endif

};

#endif /* SIMULATION_SIMULATION_CHARM_HPP */
