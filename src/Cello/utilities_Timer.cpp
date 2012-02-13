// See LICENSE_CELLO file for license and copyright information

/// @file     utilities_Timer.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2008-04-24
/// @brief    Implementation of the Timer class

#include "utilities.hpp"

//----------------------------------------------------------------------

Timer::Timer () throw()
  : time_(0),
    is_running_(false)
{
}

//----------------------------------------------------------------------

void Timer::start() throw()
{ 
  is_running_ = true;
  gettimeofday(&t1_, &tz_);
}

//----------------------------------------------------------------------

float Timer::stop() throw()
{ 
  if (is_running_) {
    gettimeofday(&t2_, &tz_);
    time_ += (t2_.tv_sec-t1_.tv_sec) + 1e-6*(t2_.tv_usec-t1_.tv_usec);
    is_running_ = false;
  }
  return time_;
}

//----------------------------------------------------------------------

void Timer::clear() throw()
{ 
  stop();
  gettimeofday(&t1_, &tz_);
  time_ = 0.0;
}

//----------------------------------------------------------------------

float Timer::value() const throw()
{
  if (is_running_) {
    gettimeofday((struct timeval *) &t2_, 
		 (struct timezone *)&tz_);
    return time_ + 
      (t2_.tv_sec-t1_.tv_sec) + 1e-6*(t2_.tv_usec-t1_.tv_usec);
  } else {
    return time_;
  }
}

//----------------------------------------------------------------------
