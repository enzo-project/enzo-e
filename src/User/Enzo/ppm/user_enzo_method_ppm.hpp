// $Id: method.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     user_enzo_method_ppm.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    Implementation of Enzo PPM hydro method

#ifndef USER_ENZO_METHOD_PPM_HPP
#define USER_ENZO_METHOD_PPM_HPP

#include <vector>
#include <string>
#include "method.hpp"


class EnzoMethodPpm : public Method {

  /// @class    EnzoMethodPpm
  /// @ingroup  Enzo
  /// @brief    Encapsulate Enzo's PPM hydro method

public: // interface

  /// Create a new Method

  EnzoMethodPpm() throw()
  {}

  /// Create Method

  EnzoMethodPpm(const EnzoMethodPpm &) throw()
  {}

  /// Delete a Method

  ~EnzoMethodPpm () throw()
  {}

  /// Apply the method

  virtual void apply() throw()
  {}


};

#endif /* USER_ENZO_METHOD_PPM_HPP */
