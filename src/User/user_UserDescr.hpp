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
    for (int i=0; i<user_method_.size(); i++) {
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
    ASSERT("UserDescr::set_user_control",
	   "",user_control_ == 0);
    user_control_ = create_user_control_(name_user_control);
    return user_control_;
  };

  /// Return the control method
  UserControl * user_control ()
  { return user_control_; };

  /// Return the timestep method
  UserTimestep * set_user_timestep (std::string name_user_timestep)
  {
    ASSERT("UserDescr::set_user_timestep",
	   "",user_timestep_ == 0);
    user_timestep_ = create_user_timestep_(name_user_timestep);
    return user_timestep_;
  };

  /// Return the timestep method
  UserTimestep * user_timestep ()
  { return user_timestep_; };

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

