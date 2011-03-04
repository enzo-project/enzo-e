// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoStopping.hpp
/// @author   James Bordner (jobordner@ucsd.edu) 
/// @date     Thu Apr  1 16:14:38 PDT 2010
/// @brief    [\ref Enzo] Implementation of Enzo's UserStopping

#ifndef ENZO_ENZO_STOPPING_HPP
#define ENZO_ENZO_STOPPING_HPP

//----------------------------------------------------------------------

class EnzoStopping : public Stopping {

  /// @class    EnzoStopping
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Stopping criteria EnzoStopping for Enzo-P

public: // interface

  /// Create and initialize a new EnzoStopping object
  EnzoStopping(EnzoBlock   * enzo,
	       int          stop_cycle,
	       double       stop_time);

  /// Return whether the simulation is complete
  virtual bool complete () throw();

private: // attributes

  /// Enzo descriptor, used for current cycle and time
  EnzoBlock * enzo_;

  /// Stop cycle
  int stop_cycle_;

  /// Stop time
  double stop_time_;

};

#endif /* ENZO_ENZO_STOPPING_HPP */
