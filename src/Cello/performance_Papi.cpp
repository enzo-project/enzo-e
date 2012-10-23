// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Papi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-12-02
/// @brief    Implementation of Papi class
///
/// Wrapper functions for PAPI

#include "cello.hpp"

#include "performance.hpp"

#ifndef CONFIG_USE_PAPI
#  define PAPI_NULL 0
#endif

//----------------------------------------------------------------------

Papi::Papi() throw()
  : is_started_(false),
    is_initialized_(false),
    event_set_(PAPI_NULL),
    num_counters_(0),
    names_(),
    values_()
{
}

//======================================================================

void Papi::init() throw()
{
#ifdef CONFIG_USE_PAPI
  // http://icl.cs.utk.edu/projects/papi/wiki/PAPIC:Low_Level

  int retval;

  retval = PAPI_library_init(PAPI_VER_CURRENT);

  if (retval != PAPI_VER_CURRENT && retval > 0) {
    WARNING("Papi::init","PAPI library version mismatch!");
  } else if (retval < 0) {
    WARNING("Papi::init","PAPI initialization error!");
  } else {
    is_initialized_ = true;
  }

  retval = PAPI_create_eventset(&event_set_);

  if (retval != PAPI_OK) {
    WARNING1("PAPI::init","PAPI_create_eventset returned %d",retval);
    is_initialized_ = false;
  } else {
    is_initialized_ = true;
  }

#endif
}

//----------------------------------------------------------------------

std::string Papi::name (int id) const throw()
{
}

//----------------------------------------------------------------------

long long Papi::value (int id) const throw()
{
}

//----------------------------------------------------------------------

int Papi::num_counters() const throw()
{
  return num_counters_;
}

//----------------------------------------------------------------------

int Papi::add_counter(int event) throw()
{
#ifdef CONFIG_USE_PAPI
  if (! is_initialized_) {

    WARNING1("PAPI::add_counter",
	     "Not adding counter %d since PAPI is not initialized",
	     event);
    return 0;

  } else {

    int retval = PAPI_add_event(event_set_, event);

    if (retval != PAPI_OK) {
      WARNING2("PAPI::add_counter","PAPI_add_event for %d returned %d",
	       name.c_str(),retval);
      return 0;
    } else {
      
      char name[40];
      PAPI_event_code_to_name(event,name);
      names_.push_back(name);

      values_.push_back(0);
      
      ++ num_counters_;

      return num_counters_-1;

    }
  }
#else
  return 0;
#endif
}

//----------------------------------------------------------------------


void Papi::start() throw()
{
#ifdef CONFIG_USE_PAPI
  if (is_started_) {
    WARNING("Papi::start",
	    "Counters already started");
  } else {
    int retval = PAPI_start(event_set_);

    if (retval != PAPI_OK) {
      WARNING1("PAPI::start()","PAPI_start() returned %d",retval);
    } else {
      is_started_ = true;
    }
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
    int retval = PAPI_start(event_set_);

    if (retval != PAPI_OK) {
      WARNING1("PAPI::start()","PAPI_start() returned %d",retval);
    } else {
      is_started_ = true;
    }
  }
#endif
}

// //----------------------------------------------------------------------

// float Papi::time_real() const throw()
// {
// #ifdef CONFIG_USE_PAPI
//   return time_real_;
// #else
//   return 0.0;
// #endif
// }

// //----------------------------------------------------------------------

// float Papi::time_proc() const throw()
// {
// #ifdef CONFIG_USE_PAPI
//   return time_proc_;
// #else
//   return 0.0;
// #endif
// }

// //----------------------------------------------------------------------

// long long Papi::flop_count() const throw()
// {
// #ifdef CONFIG_USE_PAPI
//   return flop_count_;
// #else
//   return 0;
// #endif
// }

// //----------------------------------------------------------------------

// float Papi::flop_rate() const throw()
// {
// #ifdef CONFIG_USE_PAPI
//   return flop_rate_;
// #else
//   return 0.0;
// #endif
// }

// //----------------------------------------------------------------------

// void Papi::print () const throw()
// {
// #ifdef CONFIG_USE_PAPI
//   Monitor * monitor = Monitor::instance();
//   monitor->print ("Performance","PAPI Time real   = %f",time_real());
//   monitor->print ("Performance","PAPI Time proc   = %f",time_proc());
//   monitor->print ("Performance","PAPI GFlop count = %f",flop_count()*1e-9);
//   monitor->print ("Performance","PAPI GFlop rate  = %f",flop_count()*1e-9 / time_real());
// #endif
// };
