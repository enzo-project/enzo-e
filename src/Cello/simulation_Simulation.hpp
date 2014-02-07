// See LICENSE_CELLO file for license and copyright information

/// @file     simulation_Simulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-10 16:14:57
/// @brief    Interface file for the Simulation class

#ifndef SIMULATION_SIMULATION_HPP
#define SIMULATION_SIMULATION_HPP

class Factory;
class FieldDescr;
class GroupProcess;
class Hierarchy;
class Monitor;
class Parameters;
class Performance;
class Problem;

#include "mesh.decl.h"
#include "simulation.decl.h"

class Simulation : public CBase_Simulation 
{
  /// @class    Simulation
  /// @ingroup  Simulation
  /// @brief    [\ref Simulation] Class specifying a simulation to run
  ///
  /// @detailed A Simulation object encapsulates a single simulation,
  /// and Simulation objects are replicated across processes.  Simulations
  /// are defined as groups in CHARM++.

public: // interface

  /// Simulation constructor


  Simulation
  ( const char *       parameter_file,
    int                n,
    const GroupProcess * group_process = 0
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

  /// Return the dimensionality of the simulation
  int dimension() const throw()
  { return dimension_; }

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

  /// Return the group process object
  const GroupProcess * group_process() const throw()
  { return group_process_; }

  /// Return the monitor object
  Monitor * monitor() const throw()
  { return monitor_; }

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

  /// Output Simulation information
  virtual void monitor_output();

  /// Output Performance information to stdout (root process data only)
  virtual void performance_output();

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

#ifdef CELLO_DEBUG
  FILE * fp_debug() { return fp_debug_; }
#endif


protected: // functions

  /// Initialize the Config object
  void initialize_config_ () throw();

  /// Initialize the Problem object
  void initialize_problem_ () throw();

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

protected: // attributes

#ifdef CELLO_DEBUG  
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

  /// Parallel group for the simulation
  GroupProcess * group_process_;

  /// Whether the group_process_ object was allocated inside Simulation
  bool is_group_process_new_;

  /// Dimension or rank of the simulation
  int  dimension_; 

  /// Current cycle
  int cycle_;

  /// Current time
  double time_;

  /// Current timestep
  double dt_;

  /// Current stopping criteria
  bool stop_;

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

  /// Monitor object
  Monitor * monitor_;

  /// AMR hierarchy
  Hierarchy * hierarchy_;
  
  /// Field descriptor
  FieldDescr * field_descr_;

};

#endif /* SIMULATION_SIMULATION_HPP */

