#ifndef PERFORMANCE_TIMER_HPP
#define PERFORMANCE_TIMER_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2008 James Bordner
 * Copyright (C) 2008 Laboratory for Computational Astrophysics
 * Copyright (C) 2008 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */

/** 
 *********************************************************************
 *
 * @file      timer.hpp
 * @brief     Timer class
 * @author    James Bordner
 * @date      Wed Apr 23 12:40:04 PDT 2008
 *
 * A Timer is used for timing pieces of code
 *
 * $Id$
 *
 *********************************************************************
 */

class Timer {

/** 
 *********************************************************************
 *
 * @class     Timer
 * @brief     A high resolution timer
 * @ingroup   Performance
 *
 * A Timer is used for timing pieces of code
 *
 *********************************************************************
 */

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// The accumulated time
  Scalar time_;
  /// Whether the timer is currently running
  bool is_running_;
  /// Struct values for gettimeofday()
  struct timeval t1_, t2_;
  struct timezone tz_;

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

public:

  /// Create the Timer object
  inline Timer();
  /// Start the timer
  inline void start();
  /// Stop the timer
  inline Scalar stop();
  /// Clear the timer
  inline void clear();
  /// Return the value of the timer
  inline Scalar value() const;
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
  int ret_val = gettimeofday(&t1_, &tz_);
}

//----------------------------------------------------------------------

/// Stop the timer
inline Scalar Timer::stop() 
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
inline Scalar Timer::value() const 
{
  if (is_running_) {
    gettimeofday((struct timeval *)&t2_, (struct timezone *)&tz_);
    return time_ + 
      (t2_.tv_sec-t1_.tv_sec) + 1e-6*(t2_.tv_usec-t1_.tv_usec);
  } else {
    return time_;
  }
}

#endif /* PERFORMANCE_TIMER_HPP */

