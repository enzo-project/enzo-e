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
  : is_initialized_(false),
    is_started_(false),
    event_set_(PAPI_NULL),
    num_events_(0),
    event_names_()

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

int Papi::num_events() const throw()
{
  return num_events_;
}

//----------------------------------------------------------------------

std::string Papi::event_name (int index_event) const throw()
{
  return event_names_[index_event];
}

//----------------------------------------------------------------------

int Papi::add_event(std::string event) throw()
{
#ifdef CONFIG_USE_PAPI
  if (! is_initialized_) {

    WARNING1("PAPI::add_event",
	     "Not adding event %s since PAPI is not initialized",
	     event.c_str());
    return 0;

  } else {

    int retval;
    int id;
    retval = PAPI_event_name_to_code ((char *)event.c_str(),&id);

    if (retval != PAPI_OK) {
      WARNING2("PAPI::add_event","PAPI_event_name_to_code %s returned %d",
	       event.c_str(),retval);
      return 0;
    }

    retval = PAPI_add_event(event_set_, id);

    if (retval != PAPI_OK) {
      WARNING2("PAPI::add_event","PAPI_add_event for %d returned %d",
	       id,retval);
      return 0;
    } else {
      
      event_names_.push_back(event);

      ++ num_events_;

      return num_events_-1;

    }
  }
#else
  return 0;
#endif
}

//----------------------------------------------------------------------


void Papi::start_events() throw()
{
#ifdef CONFIG_USE_PAPI
  int retval;

  if (! is_started_) {

    retval = PAPI_start(event_set_);

    if (retval != PAPI_OK) {
      WARNING1("PAPI::start_events()","PAPI_start() returned %d",retval);
    }

    is_started_ = true;

  } else {
    WARNING("Papi::start_events",
	    "Events already started");
  }
#endif
}

//----------------------------------------------------------------------

void Papi::stop_events() throw()
{
#ifdef CONFIG_USE_PAPI

  if ( is_started_) {

    long long * values = new long long [num_events_];
    int retval = PAPI_stop(event_set_,values);
    delete [] values; // values not used

    if (retval == PAPI_OK) {
      is_started_ = false;
    } else {
      WARNING1("PAPI::stop_events()","PAPI_stop() returned %d",retval);
    }

  } else {
    WARNING("Papi::stop_events",
	    "Events already stopped");
  }

#endif
}

//----------------------------------------------------------------------

int Papi::event_values (long long * values) const throw()
{
  PAPI_read(event_set_,values);
  return num_events_;
}


