// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PERFORMANCE_TIMER_HPP
#define PERFORMANCE_TIMER_HPP

/// @file     performance_Timer.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Apr 23 12:40:04 PDT 2008
/// @brief    [\ref Performance] Interface and implementation of the Timer class

class Timer {

  /// @class    Timer
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Simple class for timing code sections

public: // interface

  /// Create the Timer object
  inline Timer();

  /// Start the timer
  inline void start();

  /// Stop the timer
  inline float stop();

  /// Clear the timer
  inline void clear();

  /// Return the value of the timer
  inline float value() const;

  /// Write the timer information
  void write (FILE * fp = stdout) const throw();

private: // attributes

  /// The accumulated time
  float time_;
  /// Whether the timer is currently running
  bool is_running_;
  /// Struct values for gettimeofday()
  struct timeval t1_, t2_;
  struct timezone tz_;

};

//======================================================================

/// Create the Timer object
inline Timer::Timer () 
  : time_(0),
    is_running_(false)
{
}

//----------------------------------------------------------------------

/// Start the timer
inline void Timer::start() 
{ 
  is_running_ = true;
  gettimeofday(&t1_, &tz_);
}

//----------------------------------------------------------------------

/// Stop the timer
inline float Timer::stop() 
{ 
  if (is_running_) {
    gettimeofday(&t2_, &tz_);
    time_ += (t2_.tv_sec-t1_.tv_sec) + 1e-6*(t2_.tv_usec-t1_.tv_usec);
    is_running_ = false;
  }
  return time_;
}

//----------------------------------------------------------------------

/// Clear the timer
inline void Timer::clear() 
{ 
  stop();
  gettimeofday(&t1_, &tz_);
  time_ = 0.0;
}

//----------------------------------------------------------------------

/// Return the value of the timer
inline float Timer::value() const 
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

inline void Timer::write(FILE * fp) const throw()
{
  INCOMPLETE("Papi::write()","requires Monitor::fprint()--currently only stdout");
  Monitor * monitor = Monitor::instance();

  monitor->print ("Real time = %f\n",value());

}

#endif /* PERFORMANCE_TIMER_HPP */

