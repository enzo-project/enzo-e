// $Id: enzo_ppm.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     enzo_ppm.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    Implementation of Enzo PPM hydro method

#ifndef ENZO_PPM_HPP
#define ENZO_PPM_HPP

#include <vector>
#include <string>
#include "method.hpp"

class MethodEnzoPpm : public Method {

  /// @class    MethodEnzoPpm
  /// @ingroup  Enzo
  /// @brief    Encapsulate Enzo's PPM hydro method

public: // interface

  /// Perform any method-specific initialization

  void initialize(std::string method_name) throw();

  /// Apply the method

  void apply() throw();


};

#endif /* ENZO_PPM_HPP */
