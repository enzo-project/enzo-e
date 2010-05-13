// $Id: method_MethodEnzoPpm.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_ENZO file for license and copyright information

/// @file     method_MethodEnzoPpm.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    Implementation of Enzo PPM hydro method

#ifndef METHOD_ENZO_PPM_HPP
#define METHOD_ENZO_PPM_HPP

#include <vector>
#include <string>
#include "method.hpp"

class MethodEnzoPpm : public MethodHyperbolic {

  /// @class    MethodEnzoPpm
  /// @ingroup  Enzo
  /// @brief    Encapsulate Enzo's PPM hydro method

public: // interface

  /// Perform any method-specific initialization
  void initialize() throw();

  /// Apply the method

  void apply() throw();

private: // attributes

  bool diffusion_;
  bool flattening_;
  bool steepening_; 

};

#endif /* METHOD_ENZO_PPM_HPP */
