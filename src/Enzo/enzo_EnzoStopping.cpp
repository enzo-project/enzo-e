// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoStopping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue May 11 18:06:50 PDT 2010
/// @brief    Implementation of EnzoStopping user-dependent class member functions

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoStopping::EnzoStopping 
(
 EnzoDescr * enzo,
 int         stop_cycle,
 double      stop_time
)
  : Stopping(),
    enzo_(enzo),
    stop_cycle_(stop_cycle),
    stop_time_ (stop_time)
{
}

//----------------------------------------------------------------------

bool EnzoStopping::complete () throw()
{
  ASSERT("EnzoStopping::complete",
	 "Neither Stopping::time_stop nor Stopping::cycle_stop initialized",
	 stop_time_ != -1.0 || stop_cycle_ != -1);

  int    curr_cycle = enzo_->CycleNumber;
  double curr_time  = enzo_->Time;

  return 
    ( ! ((stop_time_  == -1.0 || curr_time  < stop_time_ ) &&
	 (stop_cycle_ == -1   || curr_cycle < stop_cycle_)));
}

