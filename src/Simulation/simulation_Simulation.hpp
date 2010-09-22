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
  UserMethod * add_user_method (std::string method_name)
  { UserMethod * method = create_user_method_(method_name);
    if (method) user_method_.push_back(method); 
    return method;
  };

  /// Return the ith user method
  UserMethod * user_method (int i)
  { return user_method_.at(i); };


  /// Set the control method
  UserControl * set_user_control (std::string name_user_control)
  {
    if (user_control_ != 0) {
      WARNING_MESSAGE("UserDescr::set_user_control",
		      "Resetting user control; deleting old one");
      delete user_control_;
      user_control_ = 0;
    }

    user_control_ = create_user_control_(name_user_control);

    return user_control_;
  };

  /// Return the control method
  UserControl * user_control ()
  { return user_control_; };

  /// Return the timestep method
  UserTimestep * set_user_timestep (std::string name_user_timestep)
  {
    if (user_timestep_ != 0) {
      WARNING_MESSAGE("UserDescr::set_user_timestep",
		      "Resetting user timestep; deleting old one");
      delete user_timestep_;
      user_timestep_ = 0;
    }

    user_timestep_ = create_user_timestep_(name_user_timestep);

    return user_timestep_;
  };


  /// Initialize user classes before a simulation
  void initialize()
  {
    ASSERT("UserDescr::initialize()","user_control_ == NULL",
	   user_control_!=0);
    ASSERT("UserDescr::initialize()","user_timestep_ == NULL",
	   user_timestep_!=0);
    ASSERT("UserDescr::initialize()","user_method_.size()==0",
	   user_method_.size()>0);
    user_control_ ->initialize(data_descr_);
    user_timestep_->initialize(data_descr_);
    for (size_t i=0; i<user_method_.size(); i++) {
      user_method_[i]->initialize(data_descr_);
    }
  }

  /// Finalize user classes after a simulation
  void finalize()
  {
    ASSERT("UserDescr::finalize()","user_control_ == NULL",
	   user_control_!=0);
    ASSERT("UserDescr::finalize()","user_timestep_ == NULL",
	   user_timestep_!=0);
    ASSERT("UserDescr::finalize()","user_method_.size()==0",
	   user_method_.size()>0);
    for (size_t i=0; i<user_method_.size(); i++) {
      user_method_[i]->finalize(data_descr_);
    }
    user_timestep_->finalize(data_descr_);
    user_control_ ->finalize(data_descr_);
  }

  /// Return the timestep method
  UserTimestep * user_timestep ()
  { return user_timestep_; };

  /// Initialize user classes before advancing all blocks one cycle
  void initialize_cycle()
  {
    INCOMPLETE_MESSAGE("UserDescr::initialize_cycle","Not implemented");
  }

  /// Finalize user classes after advancing all blocks one cycle
  void finalize_cycle()
  {
    INCOMPLETE_MESSAGE("UserDescr::finalize_cycle","Not implemented");
  }

  /// Initialize user classes before advancing a block
  void initialize_block(DataBlock * data_block)
  {
    ASSERT("UserDescr::initialize_block()","user_control_ == NULL",
	   user_control_!=0);
    ASSERT("UserDescr::initialize_block()","user_timestep_ == NULL",
	   user_timestep_!=0);
    ASSERT("UserDescr::initialize_block()","user_method_.size()==0",
	   user_method_.size()>0);
    user_control_ ->initialize_block(data_block);
    user_timestep_->initialize_block(data_block);
    for (size_t i=0; i<user_method_.size(); i++) {
      user_method_[i]->initialize_block(data_block);
    }
  }

  /// Finalize user classes after advancing a block
  void finalize_block(DataBlock * data_block)
  {
    ASSERT("UserDescr::finalize_block()","user_control_ == NULL",
	   user_control_!=0);
    ASSERT("UserDescr::finalize_block()","user_timestep_ == NULL",
	   user_timestep_!=0);
    ASSERT("UserDescr::finalize_block()","user_method_.size()==0",
	   user_method_.size()>0);
    for (size_t i=0; i<user_method_.size(); i++) {
      user_method_[i]->finalize_block(data_block);
    }
    user_timestep_->finalize_block(data_block);
    user_control_ ->finalize_block(data_block);
  }


protected: // virtual functions


  /// APPLICATION INHERITENCE OVERRIDE
  /// Read user parameters and initialize user objects
  virtual void user_initialize_() throw ()
  { 
    WARNING_MESSAGE ("Simulation::user_initialize_",
		     "No user_initialize_() implementation!");
  };


  /// APPLICATION INHERITENCE OVERRIDE
  /// Create named control method.
  virtual UserControl * create_user_control_ (std::string name_user_control)
  { 
    WARNING_MESSAGE ("Simulation::create_user_control",
		     "No create_user_control() implementation!");
    return NULL;
  };

  /// APPLICATION INHERITENCE OVERRIDE
  /// Create named timestep method.
  virtual UserTimestep * create_user_timestep_ (std::string name_user_timestep)
  { 
    WARNING_MESSAGE ("Simulation::create_user_timestep",
		     "No create_user_timestep() implementation!");
    return NULL;
  };


  /// APPLICATION INHERITENCE OVERRIDE
  /// Create named user method.
  virtual UserMethod * create_user_method_ (std::string name_user_method)
  { 
    WARNING_MESSAGE ("Simulation::create_user_method",
		     "No create_user_method() implementation!");
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

  /// User control

  UserControl *             user_control_;

  /// User time-step computation

  UserTimestep *            user_timestep_;

  /// User methods

  std::vector<UserMethod *> user_method_;

  /// Data descriptions
  DataDescr * data_descr_;

};

#endif /* SIMULATION_SIMULATION_HPP */

