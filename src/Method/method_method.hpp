// $Id: method.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef METHOD_METHOD_HPP
#define METHOD_METHOD_HPP

#include <vector>
#include <string>

enum enum_argument_type {
  enum_argument_type_unknown,
  enum_argument_type_field,
  enum_argument_type_particles
};

enum enum_access_type {
  enum_access_type_unknown,
  enum_access_type_read,
  enum_access_type_write,
  enum_access_type_read_write
};


/// @file     method.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    Declaration for the Method component

class Method {

  /// @class    Method
  /// @ingroup  Method
  /// @brief    Encapsulate external method / analysis / visualization function.

public: // interface

  /// Create a new Method

  Method() throw();

  /// Create Method

  Method(const Method &) throw();

  /// Delete a Method

  ~Method () throw();

  /// Initialize the method

  virtual void initialize(std::string method_name) throw() = 0; 


  /// Apply the method

  virtual void apply() throw() = 0; 


protected: // functions

  /// Specify a Field or Particle type, and its access type

  void add_argument(enum_argument_type argument_type,
		    std::string        argument_name,
		    enum_access_type   access_type) throw(); 

protected: // attributes

  /// List of argument types, e.g. enum_argument_type_field

  std::vector<enum_argument_type> argument_types_;

  /// List of argument names, e.g. "Density", "Velocity-X", etc.

  std::vector<std::string>        argument_names_;

  /// List of argument access types, e.g. enum_access_type_read_write

  std::vector<enum_access_type>   access_types_;

};

#endif /* METHOD_METHOD_HPP */
