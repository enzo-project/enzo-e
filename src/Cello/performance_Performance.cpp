// See LICENSE_CELLO file for license and copyright information

/// @file      performance_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2009-10-15
/// @brief     Implementation of the Performance class

#include "cello.hpp"

#include "performance.hpp"

Performance::Performance ()
  : num_counters_total_(0),
    num_counters_user_(0),
    counter_names_(),
    counter_values_(),
    index_time_real_(0),
    num_regions_(0),
    region_names_(),
    region_counters_start_(),
    region_counters_(),
    region_index_()
{

  // Add non-user counters
  new_counter("time-usec");


  num_counters_user_ = 0;

  new_region("cello");

  papi_.init();
  papi_.add_event ("PAPI_FP_OPS");

}

//----------------------------------------------------------------------

Performance::~Performance()
{
}

//----------------------------------------------------------------------

void Performance::begin() throw()
{
  for (int i=0; i<num_regions_; i++) {
    region_counters_[i].resize(num_counters_total_);
    region_counters_start_[i].resize(num_counters_total_);
  }

  index_time_real_ = num_counters_user_;

  papi_.start_events();
}

//----------------------------------------------------------------------

void Performance::end() throw()
{
  papi_.stop_events();
}

//----------------------------------------------------------------------

int Performance::new_counter(std::string counter_name)
{
  ++ num_counters_total_;
  ++ num_counters_user_;

  counter_names_.push_back(counter_name);
  counter_values_.push_back(0);

  return num_counters_user_ - 1;
}

//----------------------------------------------------------------------

long long Performance::counter(int index_counter)
{
  if (index_counter < num_counters_user_) {
    return counter_values_[index_counter];
  } else if (index_counter == index_time_real_) {
    return time_real_();
  }
}

//----------------------------------------------------------------------

void Performance::assign_counter(int index_counter, long long value)
{
  if (0 <= index_counter && index_counter < num_counters_user_) {
    counter_values_[index_counter] = value;
  } else {
    WARNING2 ("Performance::assign_counter",
	     "counter index %d out of range [0,%d)",
	     index_counter,num_counters_user_);
  }

}

//----------------------------------------------------------------------

void Performance::increment_counter(int index_counter, long long value)
{
  if (0 <= index_counter && index_counter < num_counters_user_) {
    counter_values_[index_counter] += value;
  } else {
    WARNING2 ("Performance::increment_counter",
	      "counter index %d out of range [0,%d)",
	      index_counter,num_counters_user_);
  }
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

int Performance::new_region (std::string region_name) throw()
{ 
  ++ num_regions_;

  region_names_.push_back(region_name);
  region_index_[region_name] = num_regions_ - 1;

  std::vector <long long> new_counters;
  region_counters_.push_back(new_counters);
  region_counters_start_.push_back(new_counters);

  return num_regions_ - 1;
}

//----------------------------------------------------------------------

void  Performance::start_region(int index_region) throw()
{
  for (int i=0; i<num_counters_total_; i++) {
    region_counters_start_[index_region][i] = counter(i);
    TRACE2 ("start_region %d  %lld",index_region,
	    region_counters_start_[index_region][i]);
  }
}

//----------------------------------------------------------------------

void  Performance::stop_region(int index_region) throw()
{
  for (int i=0; i<num_counters_total_; i++) {
    region_counters_[index_region][i] += 
      counter(i) - region_counters_start_[index_region][i];
    TRACE2 ("stop_region %d  %lld",index_region,
	    region_counters_[index_region][i]);
  }
}

//----------------------------------------------------------------------

void Performance::region_counters(int index_region, long long * counters) throw()
{
  for (int i=0; i<num_counters_total_; i++) {
    counters[i] = region_counters_[index_region][i];
  }
}

//======================================================================


