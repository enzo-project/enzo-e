// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Papi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-12-02
/// @brief    Implementation of Papi class
///
/// Wrapper functions for PAPI

#include "performance.hpp"

//----------------------------------------------------------------------

Papi::Papi() throw()
  : is_started_(false),
    time_real_total_(0),
    time_proc_total_(0),
    flop_count_total_(0),
    flop_rate_(0),
    time_real_(0),
    time_proc_(0),
    flop_count_(0)
      
{
}

//======================================================================

void Papi::start() throw()
{
#ifdef CONFIG_USE_PAPI
  if (is_started_) {
    WARNING("Papi::start",
	    "Counters already started");
  } else {
    is_started_ = true;
    float mflop_rate;
    // @@@@ MEMORY LEAK (4 bytes) r2026 @@@@
    PAPI_flops(&time_real_total_, 
	       &time_proc_total_, 
	       &flop_count_total_,
	       &mflop_rate);
    flop_rate_ = mflop_rate * 1e6;
  }
#endif
}

//----------------------------------------------------------------------

void Papi::stop() throw()
{
#ifdef CONFIG_USE_PAPI
  if (! is_started_) {
    WARNING("Papi::stop",
	    "Counters already stopped");
  } else {
    is_started_ = false;
    float mflop_rate;
    PAPI_flops(&time_real_, 
	       &time_proc_, 
	       &flop_count_,
	       &mflop_rate);
    flop_rate_ = mflop_rate * 1e6;

    time_real_  = time_real_  - time_real_total_;
    time_proc_  = time_proc_  - time_proc_total_;
    flop_count_ = flop_count_ - flop_count_total_;
  }
#endif
}

//----------------------------------------------------------------------

float Papi::time_real() const throw()
{
#ifdef CONFIG_USE_PAPI
  if (is_started_) {
    WARNING("Papi::time_real",
	    "Counters must be stopped");
    return 0.0;
  } else {
    return time_real_;
  }
#else
  return 0.0;
#endif
}

//----------------------------------------------------------------------

float Papi::time_proc() const throw()
{
#ifdef CONFIG_USE_PAPI
  if (is_started_) {
    WARNING("Papi::time_proc",
	    "Counters must be stopped");
    return 0.0;
  } else {
    return time_proc_;
  }
#else
  return 0.0;
#endif
}

//----------------------------------------------------------------------

long long Papi::flop_count() const throw()
{
#ifdef CONFIG_USE_PAPI
  if (is_started_) {
    WARNING("Papi::flop_count",
	    "Counters must be stopped");
    return 0;
  } else {
    return flop_count_;
  }
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

float Papi::flop_rate() const throw()
{
#ifdef CONFIG_USE_PAPI
  if (is_started_) {
    WARNING("Papi::flop_rate",
	    "Counters must be stopped");
    return 0.0;
  } else {
    return flop_rate_;
  }
#else
  return 0.0;
#endif
}

//----------------------------------------------------------------------

void Papi::print () const throw()
{
#ifdef CONFIG_USE_PAPI
  Monitor * monitor = Monitor::instance();
  monitor->print ("[Performance] PAPI Time real   = %f",time_real());
  monitor->print ("[Performance] PAPI Time proc   = %f",time_proc());
  monitor->print ("[Performance] PAPI GFlop count = %f",flop_count()*1e-9);
  monitor->print ("[Performance] PAPI GFlop rate  = %f",flop_count()*1e-9 / time_real());
#endif
};
