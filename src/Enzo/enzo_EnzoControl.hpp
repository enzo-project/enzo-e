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

  EnzoControl(Error      * error,
	      Monitor    * monitor,
	      Parameters * parameters,
	      EnzoDescr  * enzo);

  /// Initialize the EnzoControl object
  virtual void initialize (DataDescr * data_descr) throw();

  /// Finalize the EnzoControl object 
  virtual void finalize (DataDescr * data_descr) throw();

  /// Return whether the simulation is complete
  virtual bool is_done () throw();

  /// Initialize a Block at the beginning of an iteration
  void initialize_block (DataBlock * block) throw();

  /// Finalize a Block after an iteration
  void finalize_block (DataBlock * block) throw();

private:

  /// Enzo descriptor object
  EnzoDescr * enzo_;

};

#endif /* ENZO_ENZO_CONTROL_HPP */
