// See LICENSE_CELLO file for license and copyright information

/// @file     method_Initial.hpp 
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Mon Jul 13 11:11:47 PDT 2009 
/// @brief    [\ref Method] Declaration for the Initial component

#ifndef METHOD_INITIAL_HPP
#define METHOD_INITIAL_HPP

class Initial {

  /// @class    Initial
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate an initial conditions generator

public: // interface

  /// Create a new Initial
  Initial(int cycle, double time) throw()
    : cycle_(cycle), time_(time)
  {};

  /// Destructor
  virtual ~Initial() throw()
  {} ;

  /// Initial time
  double time() const throw() { return time_; }

  /// Initial cycle
  int cycle() const throw() { return cycle_; }

public: // virtual functions

  /// Enforce initial conditions on the given Block or Hierarchy
  virtual void enforce (Hierarchy * hierarchy,
			const FieldDescr * field_descr,
			Block * block = NULL) throw() = 0;

protected: // attributes

  /// Initial cycle number
  int cycle_;

  /// Initial time
  double time_;

};

#endif /* METHOD_INITIAL_HPP */
