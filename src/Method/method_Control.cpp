// $Id: method_Control.cpp 1688 2010-08-03 22:34:22Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_Control.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Jan 21 16:00:00 PST 2011
/// @brief    Implementation of the Control class

#include "cello.hpp"

#include "method.hpp"

///----------------------------------------------------------------------

bool Control::is_done (int cycle, double time) throw()
{
  ASSERT ("Control::is_done",
	  "Either stopping time or cycle must be defined",
	  cycle != -1 || time != -1);

  bool cycle_converged = (cycle != -1) ? (cycle >= cycle_stop_) : false;
  bool time_converged  = (time != -1)  ? (time  >= time_stop_)  : false;

  return (cycle_converged || time_converged);
}

///----------------------------------------------------------------------
