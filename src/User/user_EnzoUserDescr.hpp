// $Id: user_EnzoUserDescr.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef USER_ENZO_USER_DESCR_HPP
#define USER_ENZO_USER_DESCR_HPP

/// @file     user_EnzoUserDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Declaration of the EnzoUserDescr class

class EnzoUserDescr : public UserDescr {

  /// @class    EnzoUserDescr
  /// @ingroup  User
  /// @brief    Top-level description of user-implemented components

public: // interface

  /// Constructor
  EnzoUserDescr(Global * global) throw()
    : UserDescr(global),
      enzo_(new Enzo)
  {};

protected: // functions

  /// Create named control method.
  UserControl * create_user_control_ (std::string name_user_control);

  /// Create named timestep method.
  UserTimestep * create_user_timestep_ (std::string name_user_timestep);

  /// Create named user method.
  UserMethod * create_user_method_ (std::string name_user_method);

private: // attributes

  Enzo * enzo_;

};

#endif /* USER_ENZO_USER_DESCR_HPP */

