// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_Simulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-10 16:14:57
/// @brief    Interface file for the Simulation class

#ifndef SIMULATION_SIMULATION_HPP
#define SIMULATION_SIMULATION_HPP

class Factory;
class FieldDescr;
class ParticleDescr;
class ScalarDescr;
class Hierarchy;
class Monitor;
class Parameters;
class Performance;
class Problem;
class Schedule;

#include <errno.h>
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
    int                n  );

  //==================================================
  // CHARM
  //==================================================

  /// Initialize an empty Simulation
  Simulation();

  /// Initialize a migrated Simulation
  Simulation (CkMigrateMessage *m);

  //==================================================

  /// Destructor
  virtual ~Simulation();

  /// CHARM++ Pack / Unpack function
  virtual void pup (PUP::er &p);

  //----------------------------------------------------------------------
  // BLOCK INITIALIZATION WITH MsgRefine
  //----------------------------------------------------------------------

  /// Request by newly created Block to get its MsgRefine object
  virtual void p_get_msg_refine(Index index);

  /// Set MsgRefine * for a newly created Block
  void set_msg_refine (Index index, MsgRefine *);
  /// Return MsgRefine * for a newly created Block and remove from list
  MsgRefine * get_msg_refine (Index index);

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

  /// Return the scalar descriptors
  ScalarDescr * scalar_descr_long_double() throw()
  { return scalar_descr_long_double_; }
  ScalarDescr * scalar_descr_double() throw()
  { return scalar_descr_double_; }
  ScalarDescr * scalar_descr_int() throw()
  { return scalar_descr_int_; }
  ScalarDescr * scalar_descr_sync() throw()
  { return scalar_descr_sync_; }
  ScalarDescr * scalar_descr_void() throw()
  { return scalar_descr_void_; }

  /// Return the field descriptor
  FieldDescr * field_descr() const throw()
  { return field_descr_; }

  /// Return the particle descriptor
  ParticleDescr * particle_descr() const throw()
  { return particle_descr_; }

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

  /// Return true iff cycle_ changes
  bool cycle_changed() {
    bool value = false;
    if (cycle_ != cycle_watch_) {
      value = true;
      cycle_watch_ = cycle_;
    }
    return value;
  }
  
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

  /// Return the load balancing schedule
  Schedule * schedule_balance() const throw() 
  { return schedule_balance_; };

  /// Write performance information to disk (all process data)
  void performance_write();

#ifdef CONFIG_USE_PROJECTIONS  
  /// Set whether performance tracing with projections is enabled or not
  void set_projections_tracing (bool value)
  { projections_tracing_ = value; }

  bool projections_tracing() const
  { return projections_tracing_; }

  Schedule * projections_schedule_on() const
  { return projections_schedule_on_; }

  Schedule * projections_schedule_off() const
  { return projections_schedule_off_; }
#endif

public: // virtual functions

  /// Update Simulation state, including cycle, time, timestep, and
  /// stopping criteria
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
  
  // Performance

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

  /// Wait for all Hierarchy to be initialized before creating any Blocks
  void r_initialize_block_array(CkReductionMsg * msg);

  /// Wait for all local patches to be created before calling run
  void r_initialize_hierarchy(CkReductionMsg * msg);

  /// Send Config and Parameters from ip==0 to all other processes

  void send_config();
  void r_recv_config(CkReductionMsg * msg);

  //--------------------------------------------------
  // NEW OUTPUT
  //--------------------------------------------------

  void new_output_start ();

  //--------------------------------------------------
  // OLD OUTPUT
  //--------------------------------------------------

  /// Call output on Problem list of Output objects
  void output_enter ();
  void p_output_start (int index_output)
  { output_start (index_output); }
  /// Barrier between creating file(s) and writing to them
  void r_output_barrier(CkReductionMsg * msg);
  
  void output_start (int index_output);
  void output_exit();

  /// Reduce output, using p_output_write to send data to writing processes
  void s_write()
  {
    performance_->start_region(perf_output);
    write_();
    performance_->stop_region(perf_output);
  };
  void write_();

  /// Continue on to Problem::output_wait()
  void r_write(CkReductionMsg * msg);

  /// Continue on to Problem::output_wait() from checkpoint
  virtual void r_write_checkpoint();

  /// Receive data from non-writing process, write to disk, close, and
  /// proceed with next output
  void p_output_write (int n, char * buffer);

  //--------------------------------------------------
  // Compute
  //--------------------------------------------------
  
  void compute ();

  //--------------------------------------------------
  // Monitor
  //--------------------------------------------------

  void monitor_output();

   void p_monitor_performance()
   { monitor_performance(); };

  void monitor_performance();

  /// Reduction for performance data
  void r_monitor_performance_reduce (CkReductionMsg * msg);

  float timer() { return timer_.value(); }
  
  //--------------------------------------------------
  // Data
  //--------------------------------------------------

  /// Set block_array proxy on all processes
  void p_set_block_array(CProxy_Block block_array);
  
  /// Add a new Block to this local branch
  void data_insert_block(Block *) ;

  /// Remove a Block from this local branch
  void data_delete_block(Block *) ;

  /// Add a new Particle to this local branch
  void data_insert_particles(int64_t count) ;

  /// Remove a Particle from this local branch
  void data_delete_particles(int64_t count) ;

  void set_checkpoint(char * checkpoint)
  { strncpy (dir_checkpoint_,checkpoint,255);}

  //--------------------------------------------------
  // New Refresh
  //--------------------------------------------------

  /// refresh_register
  int new_register_refresh (const Refresh & refresh)
  {
    const int id_refresh = new_refresh_list_.size();
    ASSERT("Simulation::new_register_refresh()",
	   "id_refresh must be >= 0",
	   (id_refresh >= 0));
    new_refresh_list_.push_back(refresh);
    new_refresh_list_[id_refresh].set_id(id_refresh);
#ifdef DEBUG_NEW_REFRESH  
    CkPrintf ("DEBUG_NEW_REFRESH register id %d\n",id_refresh);
#endif    
    return id_refresh;
  }
  void new_refresh_set_name (int id, std::string name)
  {
    if (id >= int(new_refresh_name_.size()))
      new_refresh_name_.resize(id+1);
    new_refresh_name_[id] = name;
#ifdef DEBUG_NEW_REFRESH  
    CkPrintf ("DEBUG_NEW_REFRESH register name %d %s\n",id,name.c_str());
#endif    
  }
  
  std::string new_refresh_name (int id) const
  {
    return (0 <= id && id < int(new_refresh_name_.size())) ?
      new_refresh_name_[id] : "UNKNOWN";
  }

  /// Return the given refresh object
  Refresh & new_refresh_list (int id_refresh)
  { return new_refresh_list_[id_refresh]; }

  /// Return the number of refresh objects registered
  int new_refresh_count() const
  { return new_refresh_list_.size(); }

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

  /// Initialize the array of octrees
  void initialize_block_array_ () throw();

  /// Initialize the data object
  void initialize_data_descr_ () throw();

  /// Initialize load balancing
  void initialize_balance_ () throw();

  void deallocate_() throw();

  Schedule * create_schedule_(std::string var,
			      std::string type,
			      double start,
			      double stop,
			      double step);

  void create_checkpoint_link() {
    if (CkMyPe() == 0) {
      CkPrintf ("creating symlink %s -> %s\n",
		dir_checkpoint_,
		"Checkpoint");

      unlink ("Checkpoint");
      if (symlink(dir_checkpoint_,"Checkpoint")) {
	CkPrintf("Error: symlink(%s,\"Checkpoint\") returned %s\n",
		 dir_checkpoint_,strerror(errno));
      }
    }
  }

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

  /// Cycle at last start of performance monitoring
  int cycle_watch_;

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

  /// Schedule for projections on / off

#ifdef CONFIG_USE_PROJECTIONS
  bool projections_tracing_;
  Schedule * projections_schedule_on_;
  Schedule * projections_schedule_off_;
#endif

  /// Load balancing schedule
  Schedule * schedule_balance_;

  /// Monitor object
  Monitor * monitor_;

  /// AMR hierarchy
  Hierarchy * hierarchy_;

  /// Scalar descriptors
  ScalarDescr * scalar_descr_long_double_;
  ScalarDescr * scalar_descr_double_;
  ScalarDescr * scalar_descr_int_;
  ScalarDescr * scalar_descr_sync_;
  ScalarDescr * scalar_descr_void_;
  
  /// Field descriptor
  FieldDescr * field_descr_;

  /// Particle descriptor
  ParticleDescr * particle_descr_;

  /// Output synchronization (depreciated)
  Sync sync_output_begin_;
  Sync sync_output_write_;

  Sync sync_new_output_start_;
  Sync sync_new_output_next_;

  /// Refresh phase lists

  std::vector < Refresh >     new_refresh_list_;
  std::vector < std::string > new_refresh_name_;

  /// Saved latest checkpoint directory for creating symlink
  char dir_checkpoint_[256];

  std::map<Index,MsgRefine *> msg_refine_map_;

  /// Currently active output object
  int index_output_;
};

#endif /* SIMULATION_SIMULATION_HPP */

