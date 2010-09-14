// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef USER_USER_DESCR_HPP
#define USER_USER_DESCR_HPP

/// @file     user_UserDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Declaration of the UserDescr class

class UserDescr {

  /// @class    UserDescr
  /// @ingroup  User
  /// @brief    Top-level description of user-implemented components

public: // interface

  /// Constructor
  UserDescr(Global * global) throw()
    : user_control_(0),
      user_timestep_(0),
      user_method_(0),
      global_ (global)
  {};

  /// Destructor
  ~UserDescr() throw()
  {
    delete user_control_;
    delete user_timestep_;
    for (size_t i=0; i<user_method_.size(); i++) {
      delete user_method_[i];
    }
  }

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
      WARNING_MESSAGE("UserDescr::set_user_control","Resetting user control");
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
      WARNING_MESSAGE("UserDescr::set_user_timestep","Resetting user timestep");
    }
    user_timestep_ = create_user_timestep_(name_user_timestep);
    return user_timestep_;
  };

  /// Return the timestep method
  UserTimestep * user_timestep ()
  { return user_timestep_; };


  /// Initialize user classes before a simulation
  void initialize(DataDescr * data_descr)
  {
    ASSERT("UserDescr::initialize()","user_control_ == NULL",
	   user_control_!=0);
    ASSERT("UserDescr::initialize()","user_timestep_ == NULL",
	   user_timestep_!=0);
    ASSERT("UserDescr::initialize()","user_method_.size()==0",
	   user_method_.size()>0);
    user_control_ ->initialize(data_descr);
    user_timestep_->initialize(data_descr);
    for (size_t i=0; i<user_method_.size(); i++) {
      user_method_[i]->initialize(data_descr);
    }
  }

  /// Finalize user classes after a simulation
  void finalize(DataDescr * data_descr)
  {
    ASSERT("UserDescr::finalize()","user_control_ == NULL",
	   user_control_!=0);
    ASSERT("UserDescr::finalize()","user_timestep_ == NULL",
	   user_timestep_!=0);
    ASSERT("UserDescr::finalize()","user_method_.size()==0",
	   user_method_.size()>0);
    for (size_t i=0; i<user_method_.size(); i++) {
      user_method_[i]->finalize(data_descr);
    }
    user_timestep_->finalize(data_descr);
    user_control_ ->finalize(data_descr);
  }

  /// Initialize user classes before advancing all blocks one cycle
  void initialize_cycle(DataDescr * data_descr)
  {
    INCOMPLETE_MESSAGE("UserDescr::initialize_cycle","Not implemented");
  }

  /// Finalize user classes after advancing all blocks one cycle
  void finalize_cycle(DataDescr * data_descr)
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


protected: // functions

  /// Create named control method.
  virtual UserControl * create_user_control_ (std::string name_user_control) = 0;

  /// Create named timestep method.
  virtual UserTimestep * create_user_timestep_ (std::string name_user_timestep) = 0;

  /// Create named user method.
  virtual UserMethod * create_user_method_ (std::string name_user_method) = 0;

protected: // attributes

  UserControl *             user_control_;
  UserTimestep *            user_timestep_;
  std::vector<UserMethod *> user_method_;
  Global     *              global_;
};

#endif /* USER_USER_DESCR_HPP */

