// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Stopping.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-10-09
/// @brief    Implementation of stopping criteria Stopping

#include "problem.hpp"

//----------------------------------------------------------------------

bool Stopping::complete (int    curr_cycle,
			 double curr_time) const throw()
{
  if ((stop_cycle_ == std::numeric_limits<int>::max()) &&
      (stop_time_  == std::numeric_limits<double>::max()) &&
      (stop_seconds_  == std::numeric_limits<double>::max())) {
    ERROR("Stopping::complete",
	  "No stopping criteria specified");
  }
    
  bool stop = ( ! ((stop_time_    == -1.0 || curr_time      < stop_time_ ) &&
		   (stop_seconds_ == -1.0 || timer_.value() < stop_seconds_ ) &&
		   (stop_cycle_   == -1   || curr_cycle     < stop_cycle_)));
    
  return stop;
}
