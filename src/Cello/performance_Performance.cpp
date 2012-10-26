// See LICENSE_CELLO file for license and copyright information

/// @file      performance_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2009-10-15
/// @brief     Implementation of the Performance class
///
/// Counters

#include "cello.hpp"

#include "performance.hpp"

Performance::Performance ()
  : papi_(),
    counter_names_(),
    counter_values_(),
    num_regions_(0),
    region_names_(),
    region_counters_(),
    region_started_(),
    region_index_(),
    papi_counters_(0),
    i0_basic(0),
    i0_papi(0),
    i0_user(0),
    n_basic(0),
    n_papi(0),
    n_user(0)
{

  new_counter(counter_type_basic,"time-usec");

  //  new_region("cello");

  papi_.init();

}

//----------------------------------------------------------------------

Performance::~Performance()
{
  delete [] papi_counters_;
  papi_counters_ = 0;
}

//----------------------------------------------------------------------

void Performance::begin() throw()
{
  int n = num_counters();

  for (int i=0; i<num_regions(); i++) {
    region_counters_[i].resize(n);
    region_started_[i] = false;
  }

  papi_.start_events();
  papi_counters_ = new long long [papi_.num_events()];
}

//----------------------------------------------------------------------

void Performance::end() throw()
{
  papi_.stop_events();
}

//----------------------------------------------------------------------

int Performance::new_counter
(
 counter_type type,
 std::string  counter_name
 )
{

  counter_names_[type].push_back(counter_name);
  counter_values_.push_back(0);

  // Assumes [basic, papi, user] ordering of counters

  int id;

  if (type == counter_type_basic) {
    id = type_index_to_id (counter_type_basic,n_basic);
    ++n_basic;
    ++i0_papi;
    ++i0_user;
  } else if (type == counter_type_papi) {
    id = type_index_to_id (counter_type_papi,n_papi);
    ++n_papi;
    ++i0_user;
    papi_.add_event (counter_name);
  } else if (type == counter_type_user) {
    id = type_index_to_id (counter_type_user,n_user);
    ++n_user;
  }

  return id;
}

// //----------------------------------------------------------------------

// long long Performance::counter(int id) throw()
// {
//   refresh_counters_();
//   return counter_values_[id_to_index(id)];
// }

//----------------------------------------------------------------------

void Performance::refresh_counters_() throw()
{
  papi_.event_values(papi_counters_);

  for (int i=i0_papi; i<i0_papi+n_papi; i++) {
    counter_values_[i] = papi_counters_[i-i0_papi];
  }

  counter_values_[i0_basic] = time_real_();

}

//----------------------------------------------------------------------

void Performance::assign_counter(int id, long long value)
{
  int index = id_to_index(id);

  if (i0_user <= index && index < i0_user + n_user) {
    
    counter_values_[index] = value;

  } else {

    WARNING3 ("Performance::assign_counter",
	      "counter index %d out of range [%d,%d]",
	      index,i0_user,i0_user+n_user-1);

  }

}

//----------------------------------------------------------------------

void Performance::increment_counter(int id, long long value)
{
  int index = id_to_index(id);

  if (i0_user <= index && index < i0_user + n_user) {

    counter_values_[index] += value;

  } else {

    WARNING3 ("Performance::increment_counter",
	      "counter index %d out of range [%d,%d]",
	      index,i0_user,i0_user+n_user-1);

  }
}

//----------------------------------------------------------------------

int Performance::num_regions() const throw()
{
  return region_names_.size();
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
  int region_index = num_regions();

  region_names_.push_back(region_name);
  region_index_[region_name] = region_index;

  std::vector <long long> counters;
  region_counters_.push_back(counters);
  region_started_.push_back(false);

  return region_index;
}

//----------------------------------------------------------------------

void  Performance::start_region(int id_region) throw()
{
  // NOTE: identical to stop_region()

  TRACE1("Performance::start_region %s",region_names_[id_region].c_str());

  int index_region = id_region;

  refresh_counters_();

  if (! region_started_[index_region]) {

    region_started_[index_region] = true;

  } else {
    WARNING1 ("Performance::start_region",
	     "Region %s already started",
	     region_names_[id_region].c_str());
    return;
  }

  for (int i=0; i<num_counters(); i++) {

    region_counters_[index_region][i] = 
      counter_values_[i] - region_counters_[index_region][i];

  }
}

//----------------------------------------------------------------------

void  Performance::stop_region(int id_region) throw()
{
  // NOTE: identical to start_region()

  TRACE1("Performance::stop_region %s",region_names_[id_region].c_str());

  int index_region = id_region;

  if (region_started_[index_region]) {

    region_started_[index_region] = false;

  } else {
    WARNING1 ("Performance::stop_region",
	     "Region %s already stopped",
	     region_names_[id_region].c_str());
  }

  refresh_counters_();

  for (int i=0; i<num_counters(); i++) {

    region_counters_[index_region][i] = 
      counter_values_[i] - region_counters_[index_region][i];

  }
}

//----------------------------------------------------------------------

void Performance::region_counters(int index_region, long long * counters) throw()
{
  if (!region_started_[index_region]) {
    for (int i=0; i<num_counters(); i++) {
      counters[i] = region_counters_[index_region][i];
    }
  } else {
    refresh_counters_();
    for (int i=0; i<num_counters(); i++) {
      counters[i] = counter_values_[i] - region_counters_[index_region][i];
    }
  }
}

//----------------------------------------------------------------------

int Performance::index_to_id (int index) const throw()
{
  int id;

  if (i0_user <= index && index < i0_user + n_user) {
    id = (index - i0_user) + base_user;
  } else if (i0_basic <= index && index < i0_basic + n_basic) {
    id = (index - i0_basic) + base_basic;
  } else if (i0_papi <= index && index < i0_papi + n_papi) {
    id = (index - i0_papi) + base_papi;
  } else {
    WARNING1 ("Performance::index_to_id",
	      "counter index %d out of range",
	      index);
  }

  return id;
}

//----------------------------------------------------------------------

int Performance::type_index_to_id (counter_type type, int index) const throw()
{
  int id;

  if (type == counter_type_user) {
    id = index + base_user;
  } else if (type == counter_type_basic) {
    id = index + base_basic;
  } else if (type == counter_type_papi) {
    id = index + base_papi;
  } else {
    WARNING1 ("Performance::type_index_to_id",
	      "unknown counter_type %d",
	      type);
  }

  return id;
}

//----------------------------------------------------------------------

int Performance::id_to_index(int id) const throw()
{
  int index = 0;
  if (base_user <= id && id < base_papi) {
    index = id - base_user + i0_user;
  } else if (base_papi <= id && id < base_basic) {
    index = id - base_papi + i0_papi;
  } else if (base_basic <= id) {
    index = id - base_user + i0_basic;
  }
  return index;
}

//======================================================================

void Performance::id_to_type_index_
(int id, counter_type * type, int * index) const throw()
{
  *index = id;
  if (base_user <= id && id < base_papi) {
    *type = counter_type_user;
    *index -= base_user;
  } else if (base_papi <= id && id < base_basic) {
    *type = counter_type_papi;
    *index -= base_papi;
  } else if (base_basic <= id) {
    *type = counter_type_basic;
    *index -= base_basic;
  }
}

