// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef SIMULATION_SIMULATION_HPP
#define SIMULATION_SIMULATION_HPP

/// @file     simulation_Simulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-10 16:14:57
/// @brief    Interface file for the Simulation class
/// @todo Make Simulation same CHARM++ chare type as EnzoSimulation
/// (p.89 CHARM++ manual)
/// @note     2010-12-17: code-wiki interface review

class Boundary;
class Factory;
class FieldDescr;
class GroupProcess;
class Initial;
class Mesh;
class Method;
class Monitor;
class Output;
class Parameters;
class Performance;
class Stopping;
class Timestep;


#include PARALLEL_CHARM_INCLUDE(enzo.decl.h)

#ifdef CONFIG_USE_CHARM
class Simulation : public CBase_Simulation {
#else
class Simulation {
#endif

  /// @class    Simulation
  /// @ingroup  Simulation
  /// @brief    [\ref Simulation] Class specifying a simulation to run
  ///
  /// @detailed A Simulation object encapsulates a single simulation,
  /// and Simulation objects are replicated across processes.  Simulations
  /// are defined as groups in CHARM++.

public: // interface

  /// Constructor for CHARM++

  /// Initialize the Simulation object
  Simulation
  ( const char *   parameter_file,
#ifdef CONFIG_USE_CHARM
    int            n,
    const Factory & factory,
#else
    const Factory & factory,
    GroupProcess * group_process = 0,
#endif
    int            index = 0
    ) throw();

  //==================================================
  // CHARM
  //==================================================

#ifdef CONFIG_USE_CHARM
  /// Initialize an empty Simulation
  Simulation() {TRACE("Simulation()")};

  /// Initialize a migrated Simulation
  Simulation (CkMigrateMessage *m) {TRACE("Simulation(msg)")};

  //==================================================

  /// Request all Mesh blocks to send output to main::p_output_close()
  //  void p_output(int index, int cycle, double time) throw();

  // Output
  void p_output (int cycle, double time, double dt, bool stopping) throw();

  // Monitor, test Stopping, update Boundary and ghost zones
  void refresh () throw();

  /// default reduction callback
  void p_done (CkReductionMsg * m)
  {    printf ("done(%g)\n",*((double *)m->getData()));  delete m; }

  //--------------------------------------------------
  // Output
  //--------------------------------------------------

  /// reset output index to 0
  void output_first() throw();

  /// Process the next output object if any, else proceed with simulation
  void output_next() throw();

  /// Reduce output, using p_output_write to send data to writing processes
  void p_output_reduce() throw();

  /// Receive data from non-writing process, write to disk, close, and
  /// proceed with next output
  void p_output_write (int n, char * buffer) throw();

#endif

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  virtual ~Simulation() throw();

  //----------------------------------------------------------------------
  // ACCESSOR FUNCTIONS
  //----------------------------------------------------------------------

  /// Return the dimensionality of the Simulation
  int dimension() const throw();

  /// Return the Mesh
  Mesh * mesh() const throw();
  
  /// Return the field descriptor
  FieldDescr * field_descr() const throw();

  /// Return the performance object
  Performance * performance() const throw();

  /// Return the monitor object
  Monitor * monitor() const throw();

  /// Return the stopping object, if any
  Stopping *  stopping() const throw();
  
  /// Return the time-stepping object, if any
  Timestep * timestep() const throw();

  /// Return the initialization object, if any
  Initial *  initial() const throw();

  /// Return the boundary object, if any
  Boundary * boundary() const throw();

  /// Return the number of output objects
  size_t num_output() const throw();

  /// Return the ith output object
  Output * output(int i) const throw();

  /// Return the number of methods
  size_t num_method() const throw();

  /// Return the ith method object
  Method * method(int i) const throw();

  /// Return the factory object
  const Factory * factory () const throw();

  /// Return the current cycle number
  int cycle() const throw() {return cycle_;};

  /// Return the current time
  double time() const throw() {return time_;};

  /// Return the current dt (stored from main)
  double dt() const throw() {return dt_;};

  /// Return the currint stopping criteria (stored from main reduction)
  bool stop() const throw() {return stop_; };

  /// Return the Simulation index
  size_t index() const throw() {return index_; };

  /// Return whether to update all faces at once or axis-by-axis
  bool temp_update_all() const throw() { return temp_update_all_; };

  /// Return whether to include orthogonal ghosts in updates
  bool temp_update_full() const throw() { return temp_update_full_; };

  void update_cycle(int cycle, int time, double dt, double stop) {
    cycle_ = cycle;
    time_  = time;
    dt_    = dt;
    stop_  = stop;
  };

public: // virtual functions

  /// initialize the Simulation given a parameter file
  virtual void initialize() throw();

  /// Finalize the Simulation after running it
  virtual void finalize() throw();

  /// Run the simulation
  virtual void run() throw();

  /// Load a Simulation from disk
  virtual void read() throw();

  /// Write a Simulation state to disk
  virtual void write() const throw();

protected: // functions

  /// Initialize global simulation parameters
  void initialize_simulation_ () throw();

  /// Initialize the mesh object
  void initialize_mesh_ () throw();

  /// Initialize the data object
  void initialize_data_ () throw();


  /// Initialize the stopping object
  void initialize_stopping_ () throw();

  /// Initialize the timestep object
  void initialize_timestep_() throw();

  /// Initialize the initial conditions object
  void initialize_initial_ () throw();

  /// Initialize the boundary conditions object
  void initialize_boundary_() throw();

  /// Initialize the output objects
  void initialize_output_  () throw();

  /// Initialize the method objects
  void initialize_method_  () throw();

  /// Initialize parallelism-related parameters
  void initialize_parallel_  () throw();

  void deallocate_() throw();

protected: // abstract virtual functions

  /// Create named stopping object
  virtual Stopping * 
  create_stopping_ (std::string name) throw ();

  /// Create named timestep object
  virtual Timestep * 
  create_timestep_ (std::string name) throw ();

  /// Create named initialization object
  virtual Initial * 
  create_initial_ (std::string name) throw ();

  /// Create named boundary object
  virtual Boundary * 
  create_boundary_ (std::string name) throw ();

  /// Create named output object
  virtual Output * 
  create_output_ (std::string name) throw ();

  /// Create named method object
  virtual Method * 
  create_method_ (std::string name) throw ();

protected: // attributes

  //----------------------------------------------------------------------
  // SIMULATION PARAMETERS
  //----------------------------------------------------------------------

  /// Parameters associated with this simulation
  Parameters * parameters_;

  /// Factory for creating related families of Meshes, Patches and Blocks 
  /// [abstract factory design pattern]
  const Factory * factory_; 

  /// Parallel group for the simulation
  GroupProcess * group_process_;

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

  /// [TEMPORARY: DEBUGGING] Whether to update all faces at once or axis-by-axis
  bool temp_update_all_;

  /// [TEMPORARY: DEBUGGING] Whether to include orthogonal ghosts in updates
  bool temp_update_full_;

  //----------------------------------------------------------------------
  // SIMULATION COMPONENTS
  //----------------------------------------------------------------------

  /// Index of this simulation in an ensemble
  size_t index_;

  /// Performance object
  Performance * performance_;

  /// Monitor object
  Monitor * monitor_;

  /// AMR mesh
  Mesh * mesh_;
  
  /// Field descriptor
  FieldDescr * field_descr_;

  /// Stopping criteria
  Stopping * stopping_;

  /// Time-step computation
  Timestep * timestep_;

  /// Initial conditions object
  Initial * initial_;

  /// Boundary conditions object
  Boundary * boundary_;

  /// Output objects
  std::vector<Output *> output_list_;

  /// List of method objects
  std::vector<Method *> method_list_;


  //--------------------------------------------------
  // Output
  //--------------------------------------------------

#ifdef CONFIG_USE_CHARM

  /// Index of currently active output object
  size_t index_output_;

#endif

};

#endif /* SIMULATION_SIMULATION_HPP */

