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

  EnzoControl(Error     * error,
	      Monitor   * monitor,
	      EnzoDescr * enzo)
    : Control(error,monitor),
      enzo_(enzo)
  {};

  /// Perform any global initialization independent of specific method

  void initialize (DataDescr * data_descr) throw();

  /// Perform any global finalization independent of specific method

  void finalize (DataDescr * data_descr) throw();

  /// Perform any method-independent initialization before a block is updated

  void initialize_block (DataBlock * data_block) throw();

  /// Perform any method-independent finalization after a block is updated

  void finalize_block (DataBlock * data_block) throw();

  /// Refresh a subset of ghost zones given neighboring block face(s)

  void refresh_ghost(DataBlock * data_block, 
		     bool xm=true, bool xp=true, 
		     bool ym=true, bool yp=true, 
		     bool zm=true, bool zp=true) throw();

private:
  
  EnzoDescr * enzo_;
};

#endif /* ENZO_ENZO_CONTROL_HPP */
