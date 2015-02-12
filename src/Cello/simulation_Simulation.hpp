// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_Simulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-10 16:14:57
/// @brief    Interface file for the Simulation class

#ifndef SIMULATION_SIMULATION_HPP
#define SIMULATION_SIMULATION_HPP

class Factory;
class FieldDescr;
class Hierarchy;
class Monitor;
class Parameters;
class Performance;
class Problem;
class Schedule;

#include "mesh.decl.h"
#include "simulation.decl.h"

class Simulation : public CBase_Simulation 
{
  /// @class    Simulation
  /// @ingroup  Simulation
  /// @brief    [\ref Simulation] Class specifying a simulation to run
  ///
  /// @details A Simulation object encapsulates a single simulation,
  /// and Simulation objects are replicated across processes.  Simulations
  /// are defined as groups in CHARM++.

public: // interface

  /// Simulation constructor


  Simulation
  ( const char *       parameter_file,
    int                n
    ) throw();

  //==================================================
  // CHARM
  //==================================================

   /// Initialize an empty Simulation
   Simulation();

   /// Initialize a migrated Simulation
   Simulation (CkMigrateMessage *m);

  //==================================================

  /// Destructor
  virtual ~Simulation() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  //----------------------------------------------------------------------
  // ACCESSOR FUNCTIONS
  //----------------------------------------------------------------------

  /// Return the rank of the simulation
  int rank() const throw()
  { return rank_; }

  /// Return the Problem container object
  Problem *  problem() const throw()
  { return problem_; }

  /// Return the Hierarchy
  Hierarchy * hierarchy() const throw()
  { return hierarchy_; }
  
  /// Return the Parameters
  Parameters * parameters() const throw()
  { return parameters_; }

  /// Return the Config
  const Config * config() const throw()
  { return config_; }
  
  /// Return the field descriptor
  FieldDescr * field_descr() const throw()
  { return field_descr_; }

  /// Return the performance object associated with each cycle
  Performance * performance() throw()
  { return performance_; }

  /// Return the monitor object
  Monitor * monitor() const throw()
  { return monitor_; }

  void set_cycle(int cycle) throw()
  { cycle_ = cycle; }
  void set_time(double time) throw()
  { time_ = time; }
  void set_dt(double dt) throw()
  { dt_ = dt; }
  void set_stop(bool stop) throw()
  { stop_ = stop; }

  /// Return the current cycle number
  int cycle() const throw() 
  { return cycle_; };

  /// Return the current time
  double time() const throw() 
  { return time_; };

  /// Return the current dt (stored from main)
  double dt() const throw() 
  { return dt_; };

  /// Return the current stopping criteria (stored from main reduction)
  bool stop() const throw() 
  { return stop_; };

  /// Return the current phase of the simulation
  int phase() const throw() 
  { return phase_; };

  /// Return the current phase of the simulation
  void set_phase(int phase) const throw() 
  { phase_ = phase; };

  // /// Output Performance information to stdout (root process data only)
  // virtual void performance_output();

  /// Write performance information to disk (all process data)
  void performance_write();

public: // virtual functions

  /// Update Simulation state, including cycle, time, timestep, and stopping criteria
  virtual void update_state(int cycle, double time, double dt, double stop) ;

  /// initialize the Simulation given a parameter file
  virtual void initialize() throw();

  /// Finalize the Simulation after running it
  virtual void finalize() throw();

  /// Run the simulation
  virtual void run() throw()
  { ERROR ("Simulation::run","Implictly abstract function called"); }

  /// Return a Hierarchy factory object
  virtual const Factory * factory () const throw();
  
  int & perf_counter(int perf_region) {return perf_count_[perf_region]; }

#if defined(CELLO_DEBUG) || defined(CELLO_VERBOSE)
  FILE * fp_debug() { return fp_debug_; }
#endif

  void debug_open() {
#if defined(CELLO_DEBUG) || defined(CELLO_VERBOSE)
    char buffer[40];
    sprintf(buffer,"out.debug.%03d-%03d",CkMyPe(),cycle_);
    fp_debug_ = fopen (buffer,"w");
#endif
  }

  void debug_close() {
#if defined(CELLO_DEBUG) || defined(CELLO_VERBOSE)
    fclose(fp_debug_);
    fp_debug_ = 0;
#endif
  }

  // /// Run the simulation
  // virtual void run() throw();

  /// Add a new Block to this local branch
  void insert_block() ;

  /// Remove a Block from this local branch
  void delete_block() ;

  /// Wait for all Hierarchy to be initialized before creating any Blocks
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

protected: // functions

  /// Initialize the Config object
  void initialize_config_ () throw();

  /// Initialize the Problem object
  void initialize_problem_ () throw();

  /// Initialize the Memory object
  void initialize_memory_ () throw();

  /// Initialize global simulation parameters
  void initialize_simulation_ () throw();

  /// Initialize performance objects
  void initialize_performance_ () throw();

  /// Initialize output Monitor object
  void initialize_monitor_ () throw();

  /// Initialize the hierarchy object
  void initialize_hierarchy_ () throw();

  /// Initialize the forest of octrees
  void initialize_forest_ () throw();

  /// Initialize the data object
  void initialize_data_descr_ () throw();

  void deallocate_() throw();

  Schedule * create_schedule_(std::string var,
			      std::string type,
			      double start,
			      double stop,
			      double step);

protected: // attributes

#if defined(CELLO_DEBUG) || defined(CELLO_VERBOSE)
  FILE * fp_debug_;
#endif

  //----------------------------------------------------------------------
  // SIMULATION PARAMETERS
  //----------------------------------------------------------------------

  /// Factory, created on first call to factory()
  mutable Factory * factory_;

  /// Parameters associated with this simulation
  Parameters * parameters_;

  /// Parameter file name
  std::string parameter_file_;

  /// Rank of the simulation
  int  rank_; 

  /// Current cycle
  int cycle_;

  /// Current time
  double time_;

  /// Current timestep
  double dt_;

  /// Current stopping criteria
  bool stop_;

  /// Current phase of the cycle
  mutable int phase_;

  //----------------------------------------------------------------------
  // SIMULATION COMPONENTS
  //----------------------------------------------------------------------

  /// Configuration values, read from Parameters object
  Config * config_;

  /// Problem container object
  Problem * problem_;

  /// Simulation timer
  Timer timer_;

  /// Simulation Performance object
  Performance * performance_;

  /// Performance file name format (requires %d for process rank)
  std::string performance_name_;

  /// Processor stride for writing strict processor subset of performance data
  int performance_stride_;

  /// Counter for knowing when to call Performance start() and stop()
  int perf_count_[perf_last];

  // /// Schedule for projections on / off

  // bool projections_tracing_;
  // Schedule * projections_schedule_on_;
  // Schedule * projections_schedule_off_;

  /// Monitor object
  Monitor * monitor_;

  /// AMR hierarchy
  Hierarchy * hierarchy_;
  
  /// Field descriptor
  FieldDescr * field_descr_;

  Sync sync_output_begin_;
  Sync sync_output_write_;

};

#endif /* SIMULATION_SIMULATION_HPP */

