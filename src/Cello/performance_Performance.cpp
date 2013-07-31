// See LICENSE_CELLO file for license and copyright information

/// @file      performance_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2009-10-15
/// @brief     Implementation of the Performance class
///
/// Counters

#include "cello.hpp"

#include "performance.hpp"

Performance::Performance (Config * config)
  : papi_(),
    counter_name_(),
    counter_type_(),
    counter_values_(),
    region_name_(),
    region_counters_(),
    region_started_(),
    region_index_(),
    papi_counters_(0),
    warnings_(config ? config->performance_warnings : true)

{

  // ORDER MUST MATCH index_enum
  new_counter(counter_type_rel,"time-usec");
  new_counter(counter_type_abs,"bytes-curr");
  new_counter(counter_type_abs,"bytes-high");
  new_counter(counter_type_abs,"bytes-highest");

  papi_.init();

}

//----------------------------------------------------------------------

Performance::~Performance()
{
  delete [] papi_counters_;
  papi_counters_ = 0;
}

//----------------------------------------------------------------------

void
Performance::begin() throw()
{
  TRACE ("Performance::begin()");
  int n = num_counters();

  for (int i=0; i<num_regions(); i++) {
    region_counters_[i].resize(n);
    region_started_[i] = false;
  }

  papi_.start_events();
  papi_counters_ = new long long [papi_.num_events()];
}

//----------------------------------------------------------------------

void
Performance::end() throw()
{
  TRACE ("Performance::end()");
  papi_.stop_events();
  // stop regions?
}

//----------------------------------------------------------------------

int
Performance::new_counter ( int type, std::string  counter_name )
{

  counter_name_.push_back(counter_name);
  counter_type_.push_back(type);
  counter_values_.push_back(0);

  return counter_name_.size() - 1;
}

//----------------------------------------------------------------------

void
Performance::refresh_counters_() throw()
{
  papi_.event_values(papi_counters_);

  int ip=0;
  for (int ic=0; ic<num_counters(); ic++) {
    if (counter_type_[ic] == counter_type_papi) {
      counter_values_[ic] = papi_counters_[ip++];
    }
  }

  Memory * memory = Memory::instance();
  
  counter_values_[index_time_]          = time_real_();
  counter_values_[index_bytes_]         = memory->bytes();
  counter_values_[index_bytes_high_]    = memory->bytes_high();
  counter_values_[index_bytes_highest_] = memory->bytes_highest();

}

//----------------------------------------------------------------------

void
Performance::assign_counter(int index, long long value)
{
  if ( counter_type (index) == counter_type_user ) {
    
    counter_values_[index] = value;

  } else if (warnings_) {
      WARNING2 ("Performance::assign_counter",
		"counter %s (index %d) not a user counter",
		counter_name(index).c_str(),index);
  }

}

//----------------------------------------------------------------------

void
Performance::increment_counter(int index, long long value)
{
  if ( counter_type (index) == counter_type_user ) {

    counter_values_[index] += value;

  } else if (warnings_) {

    WARNING2 ("Performance::increment_counter",
		"counter %s (index %d) not a user counter",
		counter_name(index).c_str(),index);

  }
}

//----------------------------------------------------------------------

int
Performance::region_index (std::string name) const throw()
{
  std::map<const std::string,int>::const_iterator it;
  it=region_index_.find(name);
  if (it != region_index_.end()) {
    return it->second;
  } else {
    return -1;
  }
}

//----------------------------------------------------------------------

void
Performance::new_region (int         region_index,
			 std::string region_name) throw()
{ 
  TRACE2 ("Performance::new_region (%d %s)",region_index,region_name.c_str());
  if ((size_t)region_index >= region_name_.size()) {
    region_name_.resize(region_index+1);
  }

  region_name_[region_index] = region_name;
  region_index_[region_name]  = region_index;

  std::vector <long long> counters;
  region_counters_.push_back(counters);
  region_started_.push_back(false);
}

//----------------------------------------------------------------------

void
Performance::start_region(int id_region) throw()
{
  TRACE1 ("Performance::start_region (%d)",id_region);
  // NOTE: similar to stop_region()

  TRACE1("Performance::start_region %s",region_name_[id_region].c_str());

  int index_region = id_region;

  if (! region_started_[index_region]) {

    region_started_[index_region] = true;

  } else if (warnings_) {
    WARNING1 ("Performance::start_region",
	     "Region %s already started",
	     region_name_[id_region].c_str());
    return;
  }

  refresh_counters_();
    
  for (int i=0; i<num_counters(); i++) {

    if ( counter_type(i) == counter_type_abs ) {
      region_counters_[index_region][i] = counter_values_[i];
    } else {
      region_counters_[index_region][i] -= counter_values_[i];
    }
  }
}

//----------------------------------------------------------------------

void
Performance::stop_region(int id_region) throw()
{
  TRACE1 ("Performance::stop_region (%d)",id_region);
  // NOTE: similar to start_region()

  TRACE1("Performance::stop_region %s",region_name_[id_region].c_str());

  int index_region = id_region;

  if (region_started_[index_region]) {

    region_started_[index_region] = false;

  } else if (warnings_) {
    WARNING1 ("Performance::stop_region",
	     "Region %s already stopped",
	     region_name_[id_region].c_str());
    return;
  }

  refresh_counters_();

  for (int i=0; i<num_counters(); i++) {

    if ( counter_type(i) == counter_type_abs ) {
      region_counters_[index_region][i] = counter_values_[i];
    } else {
      region_counters_[index_region][i] += counter_values_[i];
    }

  }
}

//----------------------------------------------------------------------

bool
Performance::is_region_active(int index_region) throw()
{
  return (region_started_[index_region]);
};

//----------------------------------------------------------------------

void 
Performance::region_counters(int index_region, long long * counters) throw()
{
  if ( ! region_started_[index_region]) {
    for (int i=0; i<num_counters(); i++) {
      counters[i] = region_counters_[index_region][i];
    }
  } else {
    refresh_counters_();
    for (int i=0; i<num_counters(); i++) {
      if ( counter_type (i) == counter_type_abs ) {
	counters[i] = counter_values_[i];
      } else {
	if (is_region_active(index_region)) {
	  counters[i] = counter_values_[i] + region_counters_[index_region][i];
	} else {
	  counters[i] = counter_values_[i] - region_counters_[index_region][i];
	}
      }
    }
  }
}

//======================================================================

