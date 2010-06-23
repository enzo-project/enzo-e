// $Id: method_MethodDescr.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef METHOD_METHOD_DESCR_HPP
#define METHOD_METHOD_DESCR_HPP

/// @file     method_MethodDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Declaration of the MethodDescr class

#include "error.hpp"

class MethodDescr {

  /// @class    MethodDescr
  /// @ingroup  Method
  /// @brief    Top-level description of numerical methods

public: // interface

  /// Constructor
  MethodDescr() throw()
    : method_control_(0),
      method_timestep_(0),
      method_user_(0)
  {};

  /// Add a user method
  void add_method_user (std::string method_name)
  { MethodUser * method = create_method_user_(method_name);
    if (method) method_user_.push_back(method); 
  };

  /// Return the ith user method
  MethodUser * method_user (int i)
  { return method_user_.at(i); };


  /// Set the control method
  void set_method_control (std::string name_method_control)
  {
    ASSERT("MethodDescr::set_method_control",
	   "",method_control_ == 0);
    method_control_ = create_method_control_(name_method_control);
  };

  /// Return the control method
  MethodControl * method_control ()
  { return method_control_; };

  /// Return the timestep method
  void set_method_timestep (std::string name_method_timestep)
  {
    ASSERT("MethodDescr::set_method_timestep",
	   "",method_timestep_ == 0);
    method_timestep_ = create_method_timestep_(name_method_timestep);
  };

  /// Return the timestep method
  MethodTimestep * method_timestep ()
  { return method_timestep_; };

private: // functions

  /// Create named control method.  IMPLEMENTATION IN USER SPACE
  MethodControl * create_method_control_ (std::string name_method_control);

  /// Create named timestep method.  IMPLEMENTATION IN USER SPACE
  MethodTimestep * create_method_timestep_ (std::string name_method_timestep);

  /// Create named user method.  IMPLEMENTATION IN USER SPACE
  MethodUser * create_method_user_ (std::string name_method_user);

private: // attributes

  MethodControl *           method_control_;
  MethodTimestep *          method_timestep_;
  std::vector<MethodUser *> method_user_;

};

#endif /* METHOD_METHODDESCR_HPP */

