// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef METHOD_STOPPING_HPP
#define METHOD_STOPPING_HPP

/// @file     method_Stopping.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-05-11
/// @brief    [\ref Method] Declaration of the Stopping class

class Stopping {

  /// @class    Stopping
  /// @ingroup  Method
  /// @brief    [\ref Method] Encapsulate stopping criteria

public: // interface

  /// Constructor
  Stopping(int    stop_cycle = -1,
	   double stop_time  = -1.0) throw()
    : stop_cycle_(stop_cycle),
      stop_time_ (stop_time)
  {}

  /// Return whether the simulation is done
  virtual bool complete (int    curr_cycle,
			 double curr_time) const throw()
  {
    ASSERT("Stopping::complete",
	   "Neither Stopping::time_stop nor Stopping::cycle_stop initialized",
	   stop_time_ != -1.0 || stop_cycle_ != -1);
    return 
    ( ! ((stop_time_  == -1.0 || curr_time  < stop_time_ ) &&
	 (stop_cycle_ == -1   || curr_cycle < stop_cycle_)));
  }

  /// Return stopping cycle
  double stop_cycle () const throw()
  { return stop_cycle_; };

  /// Return stopping time
  double stop_time () const throw()
  { return stop_time_; };

protected:

  /// Stop cycle
  int stop_cycle_;

  /// Stop time
  double stop_time_;

};

#endif /* METHOD_STOPPING_HPP */

