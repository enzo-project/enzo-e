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
    num_events_(0),
    event_names_(),
    num_regions_(0),
    region_events_(),
    region_index_()

{
  insert_region_("cello");
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
    delete [] values; // values not used here

    if (retval == PAPI_OK) {
      is_started_ = true;
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

int Papi::num_regions() const throw()
{
  return num_regions_;
}

//----------------------------------------------------------------------

std::string Papi::region_name (int index_region) const throw()
{
  return region_names_[index_region];
}

//----------------------------------------------------------------------

int Papi::region_index (std::string name) const throw()
{
  return region_index_.at(name);
}

//----------------------------------------------------------------------

int Papi::add_region (std::string name_region) throw()
{ 
  insert_region_(name_region);
  return num_regions_ - 1;
}

//----------------------------------------------------------------------

void  Papi::start_region(int index_region) throw()
{  
  //  region_stack_.push(region); 
  //  std::vector<long long> new_values;
  //  values_stack_.push(new_values);
}

//----------------------------------------------------------------------

void  Papi::stop_region(int index_region) throw()
{
  // if (region != region_stack_.top() ) {
  //   WARNING2("Papi::stop_region",
  // 	     "Trying to stop region %s when active region is %s",
  // 	     region.c_str(),region_stack_.top().c_str());

  // } else {
  //   region_stack_.pop();
  //   values_stack_.pop();
  // }
}

//======================================================================

void Papi::insert_region_(std::string region) throw()
{
  region_names_.push_back(region);

  region_index_[region] = region_events_.size();

  std::vector <long long> new_events;
  new_events.resize(num_events_);
  region_events_.push_back(new_events);
}

//----------------------------------------------------------------------

const long long * Papi::values (int index_region) const throw()
{
  const long long * values = &region_events_[index_region][0];
  return values;
}


