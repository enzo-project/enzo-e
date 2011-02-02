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

  /// Initialize the Control object
  virtual void initialize (DataDescr * data_descr) throw()
  {};

  /// Finalize the Control object 
  virtual void finalize (DataDescr * data_descr) throw()
  {};

  /// Return an iterator over Blocks in a Mesh
  virtual Iterator * block_loop(Patch * patch) throw()
  { return 0; };

  /// Return whether the simulation is complete
  virtual bool is_done () throw() = 0;

protected:

  /// Error object
  Error * error_;

  /// Monitor object
  Monitor * monitor_;

};

#endif /* METHOD_CONTROL_HPP */

