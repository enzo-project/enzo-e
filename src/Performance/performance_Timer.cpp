// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Timer.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Brief description of file performance_Timer.cpp
///
/// Detailed description of file performance_Timer.cpp

#include "performance.hpp"

//----------------------------------------------------------------------

// Timer::Timer() throw ()
// {
//   INCOMPLETE("Timer::Timer");
// }

// //----------------------------------------------------------------------

// Timer::~Timer() throw ()
// {
//   INCOMPLETE("Timer::!Timer");
// }

// //----------------------------------------------------------------------

// Timer::Timer(const Timer & timer) throw ()
// /// @param     timer  Object being copied
// {
//   INCOMPLETE("Timer::Timer(Timer)");
// }

// //----------------------------------------------------------------------

// Timer & Timer::operator= (const Timer & timer) throw ()
// /// @param     timer  Source object of the assignment
// /// @return    The target assigned object
// {
//   INCOMPLETE("Timer::operator=");
//   return *this;
// }

/// Create the Timer object
Timer::Timer () throw()
  : time_(0),
    is_running_(false)
{
}

//----------------------------------------------------------------------

//Start the timer
void Timer::start() throw()
{ 
  is_running_ = true;
  gettimeofday(&t1_, &tz_);
}

//----------------------------------------------------------------------

//Stop the timer
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

//Clear the timer
void Timer::clear() throw()
{ 
  stop();
  gettimeofday(&t1_, &tz_);
  time_ = 0.0;
}

//----------------------------------------------------------------------

//Return the value of the timer
float Timer::value() const throw()
{
  if (is_running_) {
    gettimeofday((struct timeval *)&t2_, (struct timezone *)&tz_);
    return time_ + 
      (t2_.tv_sec-t1_.tv_sec) + 1e-6*(t2_.tv_usec-t1_.tv_usec);
  } else {
    return time_;
  }
}


//----------------------------------------------------------------------

void Timer::print() const throw()
{
  Monitor * monitor = Monitor::instance();

  monitor->print ("[Performance] real time = %f",value());

}

