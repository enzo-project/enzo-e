// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Papi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-12-02
/// @brief    Implementation of Papi class
///
/// Wrapper functions for PAPI

#include "cello.hpp"

#include "performance.hpp"

#ifdef CONFIG_USE_PAPI

//----------------------------------------------------------------------


Papi::Papi(bool warnings) throw()
  : is_initialized_(false),
    is_started_(false),
    event_set_(PAPI_NULL),
    num_events_(0),
    event_names_(),
    warnings_(warnings)
{
}

//======================================================================

void Papi::init() throw()
{
  int retval;

  retval = PAPI_create_eventset(&event_set_);

  if (retval != PAPI_OK) {
    if (warnings_)
      WARNING1("PAPI::init","PAPI_create_eventset returned %d",retval);
    is_initialized_ = false;
  } else {
    is_initialized_ = true;
  }

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
  if (! is_initialized_) {

    if (warnings_)
      WARNING1("PAPI::add_event",
	       "Not adding event %s since PAPI is not initialized",
	       event.c_str());
    return 0;

  } else {

    int retval;
    int id;
    retval = PAPI_event_name_to_code ((char *)event.c_str(),&id);

    if (retval != PAPI_OK) {
      if (warnings_)
	WARNING2("PAPI::add_event","PAPI_event_name_to_code %s returned %d",
		 event.c_str(),retval);
      return 0;
    }

    retval = PAPI_add_event(event_set_, id);

    if (retval != PAPI_OK) {
      if (warnings_)
	WARNING2("PAPI::add_event","PAPI_add_event for %d returned %d",
		 id,retval);
      return 0;
    } else {
      
      event_names_.push_back(event);

      ++ num_events_;

      return num_events_-1;

    }
  }
}

//----------------------------------------------------------------------


void Papi::start_events() throw()
{
  int retval;

  if (! is_started_ && num_events_ > 0 ) {

    retval = PAPI_start(event_set_);

    if (retval != PAPI_OK) {
      if (warnings_)
	WARNING1("PAPI::start_events()","PAPI_start() returned %d",retval);
    }

  } else if (is_started_) {
    if (warnings_)
      WARNING("Papi::start_events",
	      "Events already started");
  }
  is_started_ = true;
}

//----------------------------------------------------------------------

void Papi::stop_events() throw()
{

  if ( is_started_ && num_events_ > 0) {

    long long * values = new long long [num_events_];
    int retval = PAPI_stop(event_set_,values);
    delete [] values; // values not used

    if (retval != PAPI_OK) {
      if (warnings_)
	WARNING1("PAPI::stop_events()","PAPI_stop() returned %d",retval);
    }

  } else if ( ! is_started_ ) {
    if (warnings_)
      WARNING("Papi::stop_events",
	      "Events already stopped");
  }
  is_started_ = false;

}

//----------------------------------------------------------------------

int Papi::event_values (long long * values) const throw()
{
  PAPI_read(event_set_,values);
  return num_events_;
}


#endif /* CONFIG_USE_PAPI */

