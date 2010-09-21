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
  {

    // Set "default" user control and timestepping routines

  };

  /// Destructor
  ~UserDescr() throw()
  {
    delete user_control_;
    delete user_timestep_;
    for (size_t i=0; i<user_method_.size(); i++) {
      delete user_method_[i];
    }
  }


protected: // functions

  /// APPLICATION INHERITENCE OVERRIDE: Create named control method.
  virtual UserControl * create_user_control_ (std::string name_user_control) = 0;

  /// APPLICATION INHERITENCE OVERRIDE: Create named timestep method.
  virtual UserTimestep * create_user_timestep_ (std::string name_user_timestep) = 0;

  /// APPLICATION INHERITENCE OVERRIDE: Create named user method.
  virtual UserMethod * create_user_method_ (std::string name_user_method) = 0;

protected: // attributes

  UserControl *             user_control_;
  UserTimestep *            user_timestep_;
  std::vector<UserMethod *> user_method_;
  Global     *              global_;
};

#endif /* USER_USER_DESCR_HPP */

