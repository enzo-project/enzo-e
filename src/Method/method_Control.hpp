// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef METHOD_CONTROL_HPP
#define METHOD_CONTROL_HPP

/// @file     method_Control.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @todo     consolidate initialize() and initialize_block()
/// @todo     move control functionality from Method to application
/// @brief    [\ref Method] Declaration of the Control class

class Control {

  /// @class    Control
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate top-level control and description of method methods

public: // interface

  /// Constructor
  Control(Error   * error,
	  Monitor * monitor) throw()
    : error_(error),
      monitor_(monitor)
  {};

  /// Initialize the simulation
  virtual void initialize () throw() = 0;

  /// Finalize the simulation
  virtual void finalize () throw() = 0;

  /// Initialize cycle
  virtual void initialize_cycle () throw() = 0;

  /// Finalize cycle
  virtual void finalize_cycle () throw() = 0;

  /// Initialize block
  virtual void initialize_block (DataBlock * data_block) throw() = 0;

  /// Finalize block
  virtual void finalize_block (DataBlock * data_block) throw() = 0;

  /// Return whether the simulation is complete
  virtual bool is_done () throw() = 0;

protected:

  /// Error object
  Error * error_;

  /// Monitor object
  Monitor * monitor_;

};

#endif /* METHOD_CONTROL_HPP */

