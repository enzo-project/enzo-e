// See LICENSE_CELLO file for license and copyright information

/// @file      performance_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2009-10-15
/// @brief     Implementation of the Performance class

#include "cello.hpp"

#include "performance.hpp"

Performance::Performance ()
  : num_counters_(0),
    counter_names_       (),
    num_regions_(0),
    region_names_(),
    region_counters_(),
    region_index_()
{
  num_counters_ += 1; // time-real
 
  insert_region_("cello");
  papi_.start_events();
}

//----------------------------------------------------------------------

Performance::~Performance()
{
  papi_.stop_events();
}

//----------------------------------------------------------------------

void Performance::start (int index_region) throw ()
{
  timer_.start();
}

//----------------------------------------------------------------------

void Performance::stop (int index_region) throw ()
{
  timer_.stop();
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

int Performance::num_regions() const throw()
{
  return num_regions_;
}

//----------------------------------------------------------------------

std::string Performance::region_name (int index_region) const throw()
{
  return region_names_[index_region];
}

//----------------------------------------------------------------------

int Performance::region_index (std::string name) const throw()
{
  return region_index_.at(name);
}

//----------------------------------------------------------------------

int Performance::add_region (std::string name_region) throw()
{ 
  insert_region_(name_region);
  return num_regions_ - 1;
}

//----------------------------------------------------------------------

void  Performance::start_region(int index_region) throw()
{  
  //  region_stack_.push(region); 
  //  std::vector<long long> new_values;
  //  values_stack_.push(new_values);
}

//----------------------------------------------------------------------

void  Performance::stop_region(int index_region) throw()
{
  // if (region != region_stack_.top() ) {
  //   WARNING2("Performance::stop_region",
  // 	     "Trying to stop region %s when active region is %s",
  // 	     region.c_str(),region_stack_.top().c_str());

  // } else {
  //   region_stack_.pop();
  //   values_stack_.pop();
  // }
}

//======================================================================

void Performance::insert_region_(std::string region) throw()
{
  region_names_.push_back(region);

  region_index_[region] = region_counters_.size();

  std::vector <long long> new_counters;
  new_counters.resize(num_counters_);
  region_counters_.push_back(new_counters);
}

