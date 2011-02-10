// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoControl.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo's UserControl

#ifndef ENZO_ENZO_CONTROL_HPP
#define ENZO_ENZO_CONTROL_HPP


class EnzoControl : public Control {

  /// @class    EnzoControl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Method control class EnzoControl for Enzo-P

public: // interface

  /// Create a new EnzoControl

  EnzoControl(Monitor    * monitor,
	      Parameters * parameters,
	      EnzoDescr  * enzo);

  /// Initialize the simulation
  virtual void initialize () throw();

  /// Finalize the simulation
  virtual void finalize () throw();

  /// Initialize cycle
  virtual void initialize_cycle () throw();

  /// Finalize cycle
  virtual void finalize_cycle () throw();

  /// Initialize a Block at the beginning of an iteration
  virtual void initialize_block (DataBlock * block) throw();

  /// Finalize a Block after an iteration
  virtual void finalize_block (DataBlock * block) throw();

  /// Return whether the simulation is complete
  virtual bool is_done () throw();


  //--------------------------------------------------

  /// Return the Enzo descriptor
  EnzoDescr * enzo() throw() { return enzo_; };

  void read_parameters_() throw();

private: // attributes

  /// Parameters
  Parameters * parameters_;

  /// Enzo descriptor object
  EnzoDescr * enzo_;

  /// Enzo stop cycle
  int cycle_stop_;

  /// Enzo stop time
  double time_stop_;

};

#endif /* ENZO_ENZO_CONTROL_HPP */
