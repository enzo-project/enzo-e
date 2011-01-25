// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef SIMULATION_SIMULATION_HPP
#define SIMULATION_SIMULATION_HPP

/// @file     simulation_Simulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-10 16:14:57
/// @brief    Interface file for the Simulation class
/// @todo     Remove unnecessary Parameters * from function parameters
/// @note     2010-12-17: code-wiki interface review

class Simulation {

  /// @class    Simulation
  /// @ingroup  Simulation
  /// @brief    [\ref Simulation] Class specifying a simulation to run

public: // interface

  /// Initialize the Simulation object
  Simulation(Error      * error = 0,
	     Monitor    * monitor = 0,
	     Parameters * parameters = 0);

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Simulation() throw();

  /// Copy constructor
  Simulation(const Simulation & simulation) throw();

  /// Assignment operator
  Simulation & operator= (const Simulation & simulation) throw();

  //----------------------------------------------------------------------

  /// initialize the Simulation given a parameter file
  virtual void initialize(FILE * parameter_file) throw();

  /// Finalize the Simulation after running it
  void finalize() throw();

  /// Run the simulation
  void run() throw();

  /// Load a Simulation from disk
  void read() throw();

  /// Write a Simulation state to disk
  void write() throw();

  //----------------------------------------------------------------------
  // ACCESSOR FUNCTIONS
  //----------------------------------------------------------------------

  /// Return the dimensionality of the Simulation
  int dimension() throw();

  /// Return the domain extents
  void extents (double * xmin, double *xmax,
		double * ymin = 0, double *ymax = 0,
		double * zmin = 0, double *zmax = 0) throw();

  /// Return the Simulation's Error object
  Error * error() const throw();

  /// Return the Simulation's Monitor object
  Monitor * monitor() const throw();


  /// Return the Mesh
  Mesh * mesh() const throw();
  
  /// Data descriptor
  DataDescr * data_descr() const throw();

  /// Return the control method
  Control * control() const throw();
  
  /// Return the time-stepping method
  Timestep * timestep() const throw();

  /// Return the number of initialization routines
  int num_initial() const throw();

  /// Return the ith initialization routine
  Initial * initial(int i) const throw();

  /// Return the number of numerical methods
  int num_method() const throw();

  /// Return the ith method
  Method * method(int i) const throw();

protected: // functions

  /// Initialize global simulation parameters
  void initialize_simulation_ () throw();

  /// Initialize the mesh component
  void initialize_mesh_ () throw();

  /// Initialize the data component
  void initialize_data_ () throw();


  /// Initialize the control component
  void initialize_control_    () throw();

  /// Initialize the timestep component
  void initialize_timestep_   () throw();

  /// Initialize the initialization method component
  void initialize_initial_    () throw();

  /// Initialize the method components
  void initialize_method_    () throw();


protected: // abstract virtual functions

  /// Create named control method
  virtual Control * 
  create_control_ (std::string name_control) throw () = 0;

  /// Create named timestep method
  virtual Timestep * 
  create_timestep_ (std::string name_timestep) throw () = 0;

  /// Create named initialization method
  virtual Initial * 
  create_initial_ (std::string name_initial) throw () = 0;

  /// Create named numerical method
  virtual Method * 
  create_method_ (std::string name_method) throw () = 0;

protected: // attributes

  //----------------------------------------------------------------------
  // SIMULATION PARAMETERS
  //----------------------------------------------------------------------

  /// Dimension or rank of the simulation
  int  dimension_; 

  /// Lower and upper domain extents
  double extent_[6];

  //----------------------------------------------------------------------
  // SIMULATION COMPONENTS
  //----------------------------------------------------------------------

  /// Error object
  Error * error_;

  /// Monitor object
  Monitor * monitor_;

  /// Parameters associated with this simulation
  Parameters * parameters_;

  /// AMR mesh
  Mesh * mesh_;
  
  /// Data descriptor
  DataDescr * data_descr_;


  /// Method for overall control of the simulation
  Control * control_;

  /// Method for time-step computation
  Timestep * timestep_;

  /// List of initialization routines
  std::vector<Initial *> initial_list_;

  /// List of numerical methods to apply at each timestep
  std::vector<Method *> method_list_;


};

#endif /* SIMULATION_SIMULATION_HPP */

