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
    i0_basic_rel_(0),
    i0_basic_abs_(0),
    i0_papi_(0),
    i0_user_(0),
    n_basic_rel_(0),
    n_basic_abs_(0),
    n_papi_(0),
    n_user_(0)
{

  new_counter(counter_type_basic_rel,"time-usec");
  new_counter(counter_type_basic_abs,"bytes-curr");
  new_counter(counter_type_basic_abs,"bytes-high");
  new_counter(counter_type_basic_abs,"bytes-highest");

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
  // stop regions?
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

  if (type == counter_type_basic_rel) {

    id = type_index_to_id (type,n_basic_rel_);
    ++n_basic_rel_;

    ++i0_user_;    ++i0_papi_;    ++i0_basic_abs_;

  } else if (type == counter_type_basic_abs) {

    id = type_index_to_id (type,n_basic_abs_);
    ++n_basic_abs_;

    ++i0_user_;    ++i0_papi_;

  } else if (type == counter_type_papi) {

    id = type_index_to_id (type,n_papi_);
    ++n_papi_;

    ++i0_user_;

    papi_.add_event (counter_name);

  } else if (type == counter_type_user) {

    id = type_index_to_id (type,n_user_);
    ++n_user_;

  }

  TRACE4("%d %d %d %d",n_basic_abs_,n_basic_rel_,n_papi_,n_user_);

  return id;
}

//----------------------------------------------------------------------

void Performance::refresh_counters_() throw()
{
  papi_.event_values(papi_counters_);

  for (int i=i0_papi_; i<i0_papi_+n_papi_; i++) {
    counter_values_[i] = papi_counters_[i-i0_papi_];
  }

  Memory * memory = Memory::instance();
  
  counter_values_[i0_basic_rel_]   = time_real_();

  counter_values_[i0_basic_abs_]   = memory->bytes();
  counter_values_[i0_basic_abs_+1] = memory->bytes_high();
  counter_values_[i0_basic_abs_+2] = memory->bytes_highest();

}

//----------------------------------------------------------------------

void Performance::assign_counter(int id, long long value)
{
  int index = id_to_index(id);

  if ( (i0_user_ <= index) && (index < i0_user_ + n_user_ )) {
    
    counter_values_[index] = value;

  } else {

    WARNING3 ("Performance::assign_counter",
	      "counter index %d out of range [%d,%d]",
	      index,i0_user_,i0_user_+n_user_-1);

  }

}

//----------------------------------------------------------------------

void Performance::increment_counter(int id, long long value)
{
  int index = id_to_index(id);

  if ( (i0_user_ <= index) && (index < i0_user_ + n_user_) ) {

    counter_values_[index] += value;

  } else {

    WARNING3 ("Performance::increment_counter",
	      "counter index %d out of range [%d,%d]",
	      index,i0_user_,i0_user_+n_user_-1);

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
  std::map<const std::string,int>::const_iterator it;
  it=region_index_.find(name);
  if (it != region_index_.end()) {
    return it->second;
  } else {
    return -1;
  }
  //  return region_index_.at(name);
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
  // NOTE: similar to stop_region()

  TRACE1("Performance::start_region %s",region_names_[id_region].c_str());

  int index_region = id_region;

  if (! region_started_[index_region]) {

    region_started_[index_region] = true;

  } else {
    WARNING1 ("Performance::start_region",
	     "Region %s already started",
	     region_names_[id_region].c_str());
    return;
  }

  refresh_counters_();
    
  for (int i=0; i<num_counters(); i++) {

    region_counters_[index_region][i] = counter_values_[i];

  }
}

//----------------------------------------------------------------------

void  Performance::stop_region(int id_region) throw()
{
  // NOTE: similar to start_region()

  TRACE1("Performance::stop_region %s",region_names_[id_region].c_str());

  int index_region = id_region;

  if (region_started_[index_region]) {

    region_started_[index_region] = false;

  } else {
    WARNING1 ("Performance::stop_region",
	     "Region %s already stopped",
	     region_names_[id_region].c_str());
    return;
  }

  refresh_counters_();

  for (int i=0; i<num_counters(); i++) {

    if (i0_basic_abs_ <= i && i < i0_basic_abs_ + n_basic_abs_) {
      region_counters_[index_region][i] = counter_values_[i];
    } else {
      region_counters_[index_region][i] = 
	counter_values_[i] - region_counters_[index_region][i];
    }

  }
}

//----------------------------------------------------------------------

void Performance::region_counters(int index_region, long long * counters) throw()
{
  if ( ! region_started_[index_region]) {
    for (int i=0; i<num_counters(); i++) {
      counters[i] = region_counters_[index_region][i];
    }
  } else {
    refresh_counters_();
    for (int i=0; i<num_counters(); i++) {
      if (i0_basic_abs_ <= i && i < i0_basic_abs_ + n_basic_abs_) {
	counters[i] = counter_values_[i];
      } else {
	counters[i] = counter_values_[i] - region_counters_[index_region][i];
      }
    }
  }
}

//----------------------------------------------------------------------

int Performance::index_to_id (int index) const throw()
{
  int id;

  if        (i0_user_  <= index && index < i0_user_ +  n_user_) {
    id = base_user      +  (index - i0_user_);
  } else if (i0_basic_abs_ <= index && index < i0_basic_abs_ + n_basic_abs_) {
    id = base_basic_abs + (index - i0_basic_abs_);
  } else if (i0_basic_rel_ <= index && index < i0_basic_rel_ + n_basic_rel_) {
    id = base_basic_rel + (index - i0_basic_rel_);
  } else if (i0_papi_  <= index && index < i0_papi_ +  n_papi_) {
    id = base_papi      +  (index - i0_papi_);
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
  int id = index;

  if      (type == counter_type_user)  id += base_user;
  else if (type == counter_type_basic_rel) id += base_basic_rel;
  else if (type == counter_type_basic_abs) id += base_basic_abs;
  else if (type == counter_type_papi)  id += base_papi;
  else {
    WARNING1 ("Performance::type_index_to_id",
	      "unknown counter_type %d",
	      type);
  }

  return id;
}

//----------------------------------------------------------------------

int Performance::id_to_index(int id) const throw()
{
  int index = id;
  if      (base_user <= id  && id < base_user  + n_user_)
    index += (i0_user_  - base_user);
  else if (base_papi <= id  && id < base_papi  + n_papi_)  
    index += (i0_papi_  - base_papi);
  else if (base_basic_rel <= id && id < base_basic_rel + n_basic_rel_) 
    index += (i0_basic_rel_ - base_basic_rel);
  else if (base_basic_abs <= id && id < base_basic_abs + n_basic_abs_) 
    index += (i0_basic_abs_ - base_basic_abs);
  else {
    WARNING1 ("Performance::id_to_index",
	      "counter id %d out of range",
	      id);
  }
  return index;
}

//======================================================================

void Performance::id_to_type_index_
(int id, counter_type * type, int * index) const throw()
{
  (*index) = id;
  if (base_user <= id && id < base_user + n_user_) {
    (*type) = counter_type_user;
    (*index) -= base_user;
  } else if (base_papi <= id && id < base_papi + n_papi_) {
    (*type) = counter_type_papi;
    (*index) -= base_papi;
  } else if (base_basic_rel <= id && id < base_basic_rel + n_basic_rel_) {
    (*type) = counter_type_basic_rel;
    (*index) -= base_basic_rel;
  } else if (base_basic_abs <= id && id < base_basic_abs + n_basic_abs_) {
    (*type) = counter_type_basic_abs;
    (*index) -= base_basic_abs;
  } else {
    WARNING1 ("Performance::id_to_type_index_",
	      "counter id %d out of range",
	      id);
  }
    
}

