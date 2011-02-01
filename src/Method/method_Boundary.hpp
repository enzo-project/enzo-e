// $Id: method_Boundary.hpp 1896 2010-12-03 23:54:08Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Boundary.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Method] Declaration for the Boundary component

#ifndef METHOD_BOUNDARY_HPP
#define METHOD_BOUNDARY_HPP

class Boundary {

  /// @class    Boundary
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate a boundary conditions generator

public: // interface

  /// Create a new Boundary
  Boundary(Error   * error,
	   Monitor * monitor) throw()
    : error_ (error),
      monitor_ (monitor)
  {};

public: // virtual functions

  /// Perform any method-specific initialization

  virtual void initialize (DataDescr * data_descr) throw() = 0;

  /// Perform any method-specific finalization

  virtual void finalize (DataDescr * data_descr) throw() = 0;

  /// Update

  virtual void initialize_block (DataBlock * data_block) throw() = 0;

  /// Finalize PPM after advancing a block a timestep, e.g. to deallocate
  /// any dynamically-allocated variables

  virtual void finalize_block (DataBlock * data_block) throw() = 0;

  /// Return the name of the method

  virtual std::string name() const throw() = 0;

protected: // functions

  /// Specify a field or particle type and its access type

  void add_argument_(argument_enum type,
		     std::string   name,
		     access_enum   access_type,
		     DataDescr   * data_descr = 0) throw();

protected: // attributes

  /// Error
  Error * error_;

  /// Monitor
  Monitor * monitor_;

  /// List of argument types, e.g. argument_type_field
  std::vector<argument_enum> argument_types_;

  /// List of argument names, e.g. "Density", "Velocity-X", etc.
  std::vector<std::string>   argument_names_;

  /// List of argument access types, e.g. access_read_write
  std::vector<access_enum>   access_types_;

};

#endif /* METHOD_BOUNDARY_HPP */
