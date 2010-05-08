// $Id: method.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_method.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    Declaration for the Method component

#ifndef METHOD_METHOD_HPP
#define METHOD_METHOD_HPP

#include <vector>
#include <string>

#include "error.hpp"

enum argument_type {
  argument_type_unknown,
  argument_type_field,
  argument_type_particle
};

enum access_type {
  access_type_unknown,
  access_type_read,
  access_type_write,
  access_type_read_write
};


class Method {

  /// @class    Method
  /// @ingroup  Method
  /// @brief    Encapsulate external method / analysis / visualization function.

public: // interface

  /// Create a new Method
  Method() throw();

public: // virtual functions

  /// Perform any method-specific initialization

  virtual void initialize(std::string method_name) throw() = 0;

  /// Apply the method

  virtual void apply() throw() = 0; 

protected: // functions

  /// Specify a Field or Particle type and its access type

  void add_argument_(argument_type argument_type,
		     std::string   argument_name,
		     access_type   access_type) throw();

protected: // attributes

  /// Method name
  std::string method_name_;

  /// List of argument types, e.g. argument_type_field
  std::vector<argument_type> argument_types_;

  /// List of argument names, e.g. "Density", "Velocity-X", etc.
  std::vector<std::string>   argument_names_;

  /// List of argument access types, e.g. access_type_read_write
  std::vector<access_type>   access_types_;

};

#endif /* METHOD_METHOD_HPP */
