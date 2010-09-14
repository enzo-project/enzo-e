// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_USER_DESCR_HPP
#define ENZO_ENZO_USER_DESCR_HPP

/// @file     enzo_EnzoUserDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Declaration of the EnzoUserDescr class

class EnzoUserDescr : public UserDescr {

  /// @class    EnzoUserDescr
  /// @ingroup  Enzo
  /// @brief    Top-level description of user-implemented components

public: // interface

  /// Constructor
  EnzoUserDescr(Global * global) throw();

  ~EnzoUserDescr() throw()
  { delete enzo_; }

  /// Return the Enzo object created in EnzoUserDescr's constructor
  EnzoDescr * enzo() throw ()
  { return enzo_; };

protected: // functions

  /// Read parameters
  void read_parameters_() throw ();

  /// Create named control method.
  UserControl * create_user_control_ (std::string name_user_control) throw ();

  /// Create named timestep method.
  UserTimestep * create_user_timestep_ (std::string name_user_timestep) throw ();

  /// Create named user method.
  UserMethod * create_user_method_ (std::string name_user_method) throw ();

private: // attributes

  EnzoDescr * enzo_;

};

#endif /* ENZO_ENZO_USER_DESCR_HPP */

