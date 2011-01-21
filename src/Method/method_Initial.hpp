// $Id: method_MethodInitial.hpp 1896 2010-12-03 23:54:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_MethodInitial.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Method] Declaration for the MethodInitial component

#ifndef METHOD_METHOD_INITIAL_HPP
#define METHOD_METHOD_INITIAL_HPP

class MethodInitial {

  /// @class    MethodInitial
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate an initial conditions generator

public: // interface

  /// Create a new MethodInitial
  MethodInitial(Global * global) throw()
    : global_(global)
  {};

public: // virtual functions

  /// Perform any method-specific initialization

  virtual void initialize (DataDescr * data_descr) throw() = 0;

  /// Perform any method-specific finalizations steps, e.g. to
  /// deallocate any dynamically-allocated memory

  virtual void finalize (DataDescr * data_descr) throw() = 0;

  /// Initialize PPM variable that may change.  Called once per
  /// block per timestep.

  virtual void initialize_block (DataBlock * data_block) throw() = 0;

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  virtual void finalize_block (DataBlock * data_block) throw() = 0;

  /// Return the name of the method

  virtual std::string method_name() const throw() = 0;

protected: // functions

  /// Specify a field or particle type and its access type

  void add_argument_(argument_enum type,
		     std::string   name,
		     access_enum   access_type,
		     DataDescr   * data_descr = 0) throw();

protected: // attributes

  /// Global
  Global * global_;

  /// List of argument types, e.g. argument_type_field
  std::vector<argument_enum> argument_types_;

  /// List of argument names, e.g. "Density", "Velocity-X", etc.
  std::vector<std::string>   argument_names_;

  /// List of argument access types, e.g. access_read_write
  std::vector<access_enum>   access_types_;

};

#endif /* METHOD_METHOD_METHOD_HPP */
