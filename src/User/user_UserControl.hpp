// $Id: method_UserControl.hpp 1394 2010-04-22 20:52:54Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef USER_USER_CONTROL_HPP
#define USER_USER_CONTROL_HPP

/// @file     method_UserControl.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    Declaration of the UserControl class

#include <string>
#include "data.hpp"

class UserControl {

  /// @class    UserControl
  /// @ingroup  User
  /// @brief    Top-level control and description of user methods
  ///
  /// This class is not meant to be inherited from, but some functions
  /// must be implemented for user application:
  /// 
  /// 

public: // interface

  /// Constructor
  UserControl() throw()
  {};

  /// Perform any global initialization independent of specific method

  virtual void initialize (DataDescr * data_descr) throw()
  {};

  /// Perform any global finalization independent of specific method

  virtual void finalize (DataDescr * data_descr) throw()
  {};

  /// Perform any method-independent initialization before a block is updated

  virtual void initialize_block (DataBlock * data_block) throw()
  {};

  /// Perform any method-independent finalization after a block is updated

  virtual void finalize_block (DataBlock * data_block) throw()
  {};

  /// Refresh a block face's boundary / ghost zones given neighboring
  /// block face(s)

  virtual void refresh_block(DataBlock * data_block) throw()
  {};

};

#endif /* USER_USER_CONTROL_HPP */

