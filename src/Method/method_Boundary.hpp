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
  Boundary(Monitor * monitor) throw()
    : monitor_ (monitor)
  {};

public: // virtual functions

  /// Enforce boundary conditions

  virtual void initialize_block (DataBlock * data_block) throw() = 0;

  /// Perform any required clean-up

  virtual void finalize_block (DataBlock * data_block) throw() = 0;

  /// Return the name of the method

  virtual std::string name() const throw() = 0;

protected: // functions


protected: // attributes

  /// Monitor
  Monitor * monitor_;

};

#endif /* METHOD_BOUNDARY_HPP */
