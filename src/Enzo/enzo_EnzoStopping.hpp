// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoStopping.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo's UserStopping

#ifndef ENZO_ENZO_STOPPING_HPP
#define ENZO_ENZO_STOPPING_HPP


class EnzoStopping : public Stopping {

  /// @class    EnzoStopping
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Stopping criteria EnzoStopping for Enzo-P

public: // interface

  /// Create a new EnzoStopping

  EnzoStopping(Parameters * parameters,
	       EnzoDescr  * enzo);

  /// Update stopping criteria for a block
  virtual void update_block (DataBlock * block) throw();

  /// Return whether the simulation is complete
  virtual bool is_done () throw();


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

#endif /* ENZO_ENZO_STOPPING_HPP */
