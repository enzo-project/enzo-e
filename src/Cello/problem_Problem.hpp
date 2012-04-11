// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Problem.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-03-02
/// @brief    [\ref Problem] Declaration of the Problem class
///
/// This class is used as a container for classes that define the
/// problem being solved.  These classes include the following:
///
///    Boundary:    Boundary conditions
///    Initial:     Initial conditions
///    Method:      List of numerical methods
///    Output:      List of output functions
///    Refinement:  How the mesh hierarchy is to be refined
///    Stopping:    Stopping criteria
///    Timestep:    Timestepping control

#ifndef PROBLEM_PROBLEM_HPP
#define PROBLEM_PROBLEM_HPP

class Boundary;
class Factory;
class Initial;
class Method;
class Output;
class Parameters;
class Stopping;
class Timestep;

class Problem {

  /// @class    Problem
  /// @ingroup  Problem
  /// @brief    [\ref Problem] 

public: // interface

  /// Constructor
  Problem() throw();

  /// Destructor
  ~Problem() throw();

  /// Copy constructor
  Problem(const Problem & problem) throw();

  /// Assignment operator
  Problem & operator= (const Problem & problem) throw();

  /// Return the boundary object
  Boundary * boundary() const throw()  { return boundary_; }

  /// Return the ith initialization object
  Initial *  initial(int i = -1) const throw()
  {
    if (i == -1) i = index_initial_;
    return (i < (int)initial_list_.size()) ? initial_list_[i] : NULL; 
  }

#ifdef CONFIG_USE_CHARM
  /// reset initial index to 0 (not needed, but mirrors initial_output() )
  void initial_first() throw();

  /// Process the next initial object if any, else proceed with simulation
  void initial_next(Simulation * simulation) throw();

#endif


  /// Return the ith method object
  Method * method(size_t i) const throw() 
  { return (i < method_list_.size()) ? method_list_[i] : NULL; }

  /// Return the ith output object
  Output * output(int i = -1) const throw()
  { 
    if (i == -1) i = index_output_;
    return (i < (int)output_list_.size()) ? output_list_[i] : NULL; 
  }

#ifdef CONFIG_USE_CHARM
  /// reset output index to 0
  void output_first() throw();

  /// Process the next output object if any, else proceed with simulation
  void output_next(Simulation * simulation) throw();

  /// Reduce output, using p_output_write to send data to writing processes
  void output_reduce(Simulation * simulation) throw();
  
  /// Receive data from non-writing process, write to disk, close, and
  /// proceed with next output
  void output_write (Simulation * simulation, int n, char * buffer) throw();

#endif

  /// Return the stopping object
  Stopping *  stopping() const throw() { return stopping_; }

  /// Return the timestep control object
  Timestep * timestep() const throw()  { return timestep_; }

  /// Initialize the boundary conditions object
  void initialize_boundary(Parameters * parameters) throw();

  /// Initialize the initial conditions object
  void initialize_initial(Parameters * parameters,
			  GroupProcess * group_process) throw();

  /// Initialize the method objects
  void initialize_method(Parameters * parameters) throw();

  /// Initialize the output objects
  void initialize_output(Parameters * parameters,
			 FieldDescr * field_descr,
			 GroupProcess * group_process,
			 const Factory * factory) throw();

  /// Initialize the stopping object
  void initialize_stopping(Parameters * parameters) throw();

  /// Initialize the timestep object
  void initialize_timestep(Parameters * parameters) throw();


protected: // functions

  /// Deallocate components
  void deallocate_() throw();

  /// Create named boundary object
  virtual Boundary * create_boundary_
  (std::string name, Parameters * parameters) throw ();

  /// Create named initialization object
  virtual Initial *  create_initial_ 
  (std::string name, Parameters * parameters,
   GroupProcess * = 0) throw ();

  /// Create named method object
  virtual Method *   create_method_
  (std::string name, Parameters * parameters) throw ();

  /// Create named output object
  virtual Output *   create_output_  
  (std::string name, Parameters * parameters,
   GroupProcess *, const Factory * ) throw ();

  /// Create named stopping object
  virtual Stopping * create_stopping_ 
  (std::string name, Parameters * parameters) throw ();

  /// Create named timestep object
  virtual Timestep * create_timestep_ 
  (std::string name, Parameters * parameters) throw ();

private: // attributes

  /// Boundary conditions object
  Boundary * boundary_;

  /// Initial conditions object
  std::vector<Initial *> initial_list_;

  /// Stopping criteria
  Stopping * stopping_;

  /// Time-step computation
  Timestep * timestep_;

  /// List of method objects
  std::vector<Method *> method_list_;

  /// Output objects
  std::vector<Output *> output_list_;

  /// Index of currently active Initial object
  size_t index_initial_;

  /// Index of currently active Output object
  size_t index_output_;

};

#endif /* PROBLEM_PROBLEM_HPP */

