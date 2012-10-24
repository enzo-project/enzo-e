// See LICENSE_CELLO file for license and copyright information

/// @file      performance_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2009-10-15
/// @brief     Implementation of the Performance class

#include "cello.hpp"

#include "performance.hpp"

Performance::Performance ()
  : counters_(),
    counter_names_       (),
    region_names_        (),
    current_region_      (0)
{
  counters_.push_back(new PerfCounters(num_counters));
}

//----------------------------------------------------------------------

Performance::~Performance()
{
  deallocate_();
}

//----------------------------------------------------------------------

void Performance::start (int index_region) throw ()
{
  timer_.start();
  papi_.start_region(index_region);
}

//----------------------------------------------------------------------

void Performance::stop (int index_region) throw ()
{
  timer_.stop();
  papi_.stop_region(index_region);
}

//----------------------------------------------------------------------

unsigned Performance::new_region(std::string region_name)
{
  region_names_.push_back(region_name);
  return region_names_.size()-1;
}

//----------------------------------------------------------------------

int Performance::region(unsigned id_region)
{
  INCOMPLETE("Performance::region");
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_region(unsigned id_region)
{
  INCOMPLETE("Performance::set_region");
}

//----------------------------------------------------------------------

void Performance::start_region(unsigned region_id)
{
  INCOMPLETE("Performance::start_region");
}

//----------------------------------------------------------------------

void Performance::stop_region(unsigned region_id)
{
  INCOMPLETE("Performance::stop_region");
}

//----------------------------------------------------------------------

unsigned Performance::new_counter(std::string counter_name)
{
  counter_names_.push_back(counter_name);
  return counter_names_.size()-1;
}

//----------------------------------------------------------------------

type_counter Performance::counter(unsigned id_counter)
{
  INCOMPLETE("Performance::counter");
  return 0;
}

//----------------------------------------------------------------------

void Performance::set_counter(unsigned          id_counter,
			      type_counter value)
{
  INCOMPLETE("Performance::set_counter");
}

//----------------------------------------------------------------------

void Performance::increment_counter(unsigned          id_counter,
				    type_counter value)
{
  INCOMPLETE("Performance::increment_counter");
}

//----------------------------------------------------------------------

void Performance::flush()
{
  INCOMPLETE("Performance::flush");
}


//======================================================================
void Performance::deallocate_() throw()
{
  for (unsigned i=0; i<counters_.size(); i++) {
    delete counters_.at(i);
  }
}
