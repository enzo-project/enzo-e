// $Id: method.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Method.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    Declaration for the Method component

#ifndef METHOD_METHOD_HPP
#define METHOD_METHOD_HPP

#include <vector>
#include <string>

#include "error.hpp"
#include "data.hpp"

enum argument_type {
  argument_unknown,
  argument_field,
  argument_particle
};

enum access_type {
  access_unknown,
  access_read,
  access_write,
  access_read_write
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

  virtual void initialize_method(DataDescr * data_descr) throw() = 0;

  /// Perform any method-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  virtual void finalize_method(DataDescr * data_descr) throw() = 0;

  /// Initialize PPM variable that may change.  Called once per
  /// block per timestep.

  virtual void initialize_block (DataBlock * data_block) throw() = 0;

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  virtual void finalize_block (DataBlock * data_block) throw() = 0;

  /// Apply the method to advance a block one timestep 

  virtual void advance_block( DataBlock * data_block,
			      double t, double dt ) throw() = 0; 

  /// Refresh a block face's boundary / ghost zones given neighboring
  /// block face(s)

  virtual void refresh_face() throw() = 0;

protected: // functions

  /// Specify a field or particle type and its access type

  void add_argument_(argument_type type,
		     std::string   name,
		     access_type   access_type) throw();

protected: // attributes

  /// Method name
  std::string method_name_;

  /// List of argument types, e.g. argument_type_field
  std::vector<argument_type> argument_types_;

  /// List of argument names, e.g. "Density", "Velocity-X", etc.
  std::vector<std::string>   argument_names_;

  /// List of argument access types, e.g. access_read_write
  std::vector<access_type>   access_types_;

};

#endif /* METHOD_METHOD_HPP */
