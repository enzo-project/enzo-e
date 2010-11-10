// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef SIMULATION_SIMULATION_HPP
#define SIMULATION_SIMULATION_HPP

/// @file     simulation_Simulation.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2009-11-10 16:14:57
/// @brief    Interface file for the Simulation class

class Simulation {

  /// @class    Simulation
  /// @ingroup  Simulation
  /// @brief    Class specifying a simulation to run

public: // interface

  /// Initialize the Simulation object
  Simulation(Global * global);

  //----------------------------------------------------------------------
  // Big Three
  //----------------------------------------------------------------------

  /// Destructor
  ~Simulation() throw();

  /// Copy constructor
  Simulation(const Simulation & classname) throw();

  /// Assignment operator
  Simulation & operator= (const Simulation & classname) throw();

  /// Advance the simulation a specified amount

  void advance(double stop_time, int stop_cycle)  throw();

  /// Return the dimension 1 <= d <= 3, of the simulation
  int dimension() const throw(); 

  /// Return extents.  Assumes domain_lower[] and domain_upper[] are allocated to at least dimension()
  int domain (int domain_lower[], int domain_upper[]) throw();

  /// Return the Global object
  Global * global ()  throw() 
  { return global_; };

  /// Return the Mesh object
  Mesh * mesh ()  throw() 
  { return mesh_; };

  /// Return the Data code descriptor
  DataDescr * data_descr () throw() 
  { return data_descr_; };

  /// Add a user method
  MethodHyperbolic * add_method_hyperbolic (std::string method_name)
  { MethodHyperbolic * method = create_method_hyperbolic_(method_name);
    if (method) method_hyperbolic_.push_back(method); 
    return method;
  };

  /// Return the ith method hyperbolic
  MethodHyperbolic * method_hyperbolic (int i)
  { return method_hyperbolic_.at(i); };


  /// Set the control method
  MethodControl * set_method_control (std::string name_method_control)
  {
    if (method_control_ != 0) {
      WARNING_MESSAGE("MethodDescr::set_method_control",
		      "Resetting method control; deleting old one");
      delete method_control_;
      method_control_ = 0;
    }

    method_control_ = create_method_control_(name_method_control);

    return method_control_;
  };

  /// Return the control method
  MethodControl * method_control ()
  { return method_control_; };

  /// Return the timestep method
  MethodTimestep * set_method_timestep (std::string name_method_timestep)
  {
    if (method_timestep_ != 0) {
      WARNING_MESSAGE("MethodDescr::set_method_timestep",
		      "Resetting method timestep; deleting old one");
      delete method_timestep_;
      method_timestep_ = 0;
    }

    method_timestep_ = create_method_timestep_(name_method_timestep);

    return method_timestep_;
  };


  /// Initialize method classes before a simulation
  void initialize()
  {
    ASSERT("MethodDescr::initialize()","method_control_ == NULL",
	   method_control_!=0);
    ASSERT("MethodDescr::initialize()","method_timestep_ == NULL",
	   method_timestep_!=0);
    ASSERT("MethodDescr::initialize()","method_hyperbolic_.size()==0",
	   method_hyperbolic_.size()>0);
    method_control_ ->initialize(data_descr_);
    method_timestep_->initialize(data_descr_);
    for (size_t i=0; i<method_hyperbolic_.size(); i++) {
      method_hyperbolic_[i]->initialize(data_descr_);
    }
  }

  /// Finalize method classes after a simulation
  void finalize()
  {
    ASSERT("MethodDescr::finalize()","method_control_ == NULL",
	   method_control_!=0);
    ASSERT("MethodDescr::finalize()","method_timestep_ == NULL",
	   method_timestep_!=0);
    ASSERT("MethodDescr::finalize()","method_hyperbolic_.size()==0",
	   method_hyperbolic_.size()>0);
    for (size_t i=0; i<method_hyperbolic_.size(); i++) {
      method_hyperbolic_[i]->finalize(data_descr_);
    }
    method_timestep_->finalize(data_descr_);
    method_control_ ->finalize(data_descr_);
  }

  /// Return the timestep method
  MethodTimestep * method_timestep ()
  { return method_timestep_; };

  /// Initialize method classes before advancing all blocks one cycle
  void initialize_cycle()
  {
    INCOMPLETE_MESSAGE("MethodDescr::initialize_cycle","Not implemented");
  }

  /// Finalize method classes after advancing all blocks one cycle
  void finalize_cycle()
  {
    INCOMPLETE_MESSAGE("MethodDescr::finalize_cycle","Not implemented");
  }

  /// Initialize method classes before advancing a block
  void initialize_block(DataBlock * data_block)
  {
    ASSERT("MethodDescr::initialize_block()","method_control_ == NULL",
	   method_control_!=0);
    ASSERT("MethodDescr::initialize_block()","method_timestep_ == NULL",
	   method_timestep_!=0);
    ASSERT("MethodDescr::initialize_block()","method_hyperbolic_.size()==0",
	   method_hyperbolic_.size()>0);
    method_control_ ->initialize_block(data_block);
    method_timestep_->initialize_block(data_block);
    for (size_t i=0; i<method_hyperbolic_.size(); i++) {
      method_hyperbolic_[i]->initialize_block(data_block);
    }
  }

  /// Finalize method classes after advancing a block
  void finalize_block(DataBlock * data_block)
  {
    ASSERT("MethodDescr::finalize_block()","method_control_ == NULL",
	   method_control_!=0);
    ASSERT("MethodDescr::finalize_block()","method_timestep_ == NULL",
	   method_timestep_!=0);
    ASSERT("MethodDescr::finalize_block()","method_hyperbolic_.size()==0",
	   method_hyperbolic_.size()>0);
    for (size_t i=0; i<method_hyperbolic_.size(); i++) {
      method_hyperbolic_[i]->finalize_block(data_block);
    }
    method_timestep_->finalize_block(data_block);
    method_control_ ->finalize_block(data_block);
  }


protected: // virtual functions


  /// APPLICATION INHERITENCE OVERRIDE
  /// Read method parameters and initialize method objects
  virtual void method_initialize_() throw ()
  { 
    WARNING_MESSAGE ("Simulation::method_initialize_",
		     "No method_initialize_() implementation!");
  };


  /// APPLICATION INHERITENCE OVERRIDE
  /// Create named control method.
  virtual MethodControl * create_method_control_ (std::string name_method_control)
  { 
    WARNING_MESSAGE ("Simulation::create_method_control",
		     "No create_method_control() implementation!");
    return NULL;
  };

  /// APPLICATION INHERITENCE OVERRIDE
  /// Create named timestep method.
  virtual MethodTimestep * create_method_timestep_ (std::string name_method_timestep)
  { 
    WARNING_MESSAGE ("Simulation::create_method_timestep",
		     "No create_method_timestep() implementation!");
    return NULL;
  };


  /// APPLICATION INHERITENCE OVERRIDE
  /// Create named method method.
  virtual MethodHyperbolic* create_method_hyperbolic_ (std::string name_method_hyperbolic)
  { 
    WARNING_MESSAGE ("Simulation::create_method_hyperbolic",
		     "No create_method_hyperbolic() implementation!");
    return NULL;
  };

private: // functions

  /// Initialize the mesh

  void initialize_mesh_() throw ();

protected: // attributes

  /// Dimension or rank of the simulation

  int dimension_;

  /// Domain extents

  double domain_lower_[3];
  double domain_upper_[3];

  /// "global" data, including parameters, monitor, error, parallel, etc.

  Global * global_;

  /// AMR mesh
  Mesh * mesh_;

  /// Method control

  MethodControl *             method_control_;

  /// Method time-step computation

  MethodTimestep *            method_timestep_;

  /// Method methods

  std::vector<MethodHyperbolic *> method_hyperbolic_;

  /// Data descriptions
  DataDescr * data_descr_;

};

#endif /* SIMULATION_SIMULATION_HPP */

