// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef METHOD_HPP
#define METHOD_HPP

/// @file     method.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @todo     Split file into separate include and class definition files
/// @brief    Include file for the Method component

#include "cello.h"

class Method {

  /// @class    Method
  /// @ingroup  Method
  /// @brief    Encapsulate external method / analysis / visualization function.

public: // interface

  /// Create a new Method
  Method() throw();

  /// Create Method
  Method(const Method &) throw();

  /// Initialize the method
  void initialize() throw(); 

  /// Specify a Field or Particle type read or modified 
  void add_argument(std::string) throw(); 

  /// Apply the method
  void apply() throw(); 

};

#include "method_ppm.hpp"

#endif /* METHOD_HPP */
