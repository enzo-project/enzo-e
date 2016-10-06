// See LICENSE_CELLO file for license and copyright information

/// @file      performance_Performance.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      2009-10-15
/// @brief     Implementation of the Performance class
///
/// Counters

#include "cello.hpp"

#include "performance.hpp"

// #define TRACE_PERFORMANCE

static long long time_start[CONFIG_NODE_SIZE] = { 0 };

Performance::Performance (Config * config)
  :
#ifdef CONFIG_USE_PAPI  
  papi_(config ? config->performance_warnings : false),
#endif
  counter_name_(),
  counter_type_(),
  counter_values_(),
  region_name_(),
  region_counters_(),
  region_started_(),
  region_index_(),
  region_in_charm_(),
#ifdef CONFIG_USE_PAPI  
  papi_counters_(0),
#endif
  warnings_(config ? config->performance_warnings : false),
  index_region_current_(perf_unknown)
{

  const int in = cello::index_static();

  time_start[in] = time_real_();

  // ORDER MUST MATCH index_enum
  new_counter(counter_type_rel,"time-usec");
  new_counter(counter_type_abs,"bytes-curr");
  new_counter(counter_type_abs,"bytes-high");
  new_counter(counter_type_abs,"bytes-highest");

#ifdef CONFIG_USE_PAPI  
  papi_.init();
#endif  

}

//----------------------------------------------------------------------

Performance::~Performance()
{
#ifdef CONFIG_USE_PAPI  
  delete [] papi_counters_;
  papi_counters_ = NULL;
#endif  
}

//----------------------------------------------------------------------

void
Performance::begin() throw()
{
#ifdef TRACE_PERFORMANCE
  CkPrintf ("%d TRACE_PERFORMANCE Performance::begin()\n",CkMyPe());
#endif
  int n = num_counters();

  for (int i=0; i<num_regions(); i++) {
    region_counters_[i].resize(n);
    region_started_[i] = false;
  }

#ifdef CONFIG_USE_PAPI  
  papi_counters_ = new long long [papi_.num_events()];
  papi_.start_events();
#endif  
}

//----------------------------------------------------------------------

void
Performance::end() throw()
{
#ifdef TRACE_PERFORMANCE
  CkPrintf ("%d TRACE_PERFORMANCE Performance::end()\n",CkMyPe());
#endif
}

//----------------------------------------------------------------------

int
Performance::new_counter ( int type, std::string  counter_name )
{

  counter_name_.push_back(counter_name);
  counter_type_.push_back(type);
  counter_values_.push_back(0);

#ifdef CONFIG_USE_PAPI  
  if (type == counter_type_papi) {
    papi_.add_event(counter_name);
  }
#endif

  return counter_name_.size() - 1;
}

//----------------------------------------------------------------------

void
Performance::refresh_counters_() throw()
{
#ifdef CONFIG_USE_PAPI  
  papi_.event_values(papi_counters_);

  int ip=0;
  for (int ic=0; ic<num_counters(); ic++) {
    if (counter_type_[ic] == counter_type_papi) {
      counter_values_[ic] = papi_counters_[ip++];
    }
  }
#endif

  Memory * memory = Memory::instance();

  const int in = cello::index_static();

  counter_values_[index_time_]          = time_real_()-time_start[in];
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
			 std::string region_name,
			 bool        in_charm) throw()
{ 
#ifdef TRACE_PERFORMANCE
  CkPrintf ("%d TRACE_PERFORMANCE Performance::new_region (%d %s)\n",CkMyPe(),
	    region_index,region_name.c_str());
#endif
  if ((size_t)region_index >= region_name_.size()) {
    region_name_.resize(region_index+1);
    region_in_charm_.resize(region_index+1);
  }

  region_name_[region_index]    = region_name;
  region_index_[region_name]    = region_index;
  region_in_charm_[region_index] = in_charm;

  std::vector <long long> counters;
  region_counters_.push_back(counters);
  region_started_.push_back(false);
}

//----------------------------------------------------------------------

void
Performance::start_region(int id_region, std::string file, int line) throw()
{
#ifdef TRACE_PERFORMANCE
  CkPrintf ("%d TRACE_PERFORMANCE Performance::start_region (%d,%s) %s:%d\n",CkMyPe(),
	    id_region,region_name_[id_region].c_str(),file.c_str(),line);
#endif

  if (region_in_charm_[index_region_current_]) {
    stop_region(index_region_current_);
  }
  
  int index_region = id_region;

  index_region_current_ = index_region;

  if (! region_started_[index_region]) {

    region_started_[index_region] = true;

  } else if (warnings_) {
    if (file == "") {
      WARNING1 ("Performance::start_region",
		"Region %s already started",
		region_name_[id_region].c_str());
    // } else {
    //   WARNING3 ("Performance::start_region",
    // 		"Region %s already started %s %d",
    // 		region_name_[id_region].c_str(),
    // 		file.c_str(),line);
    }
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
Performance::stop_region(int id_region, std::string file, int line) throw()
{

#ifdef TRACE_PERFORMANCE
  CkPrintf ("%d TRACE_PERFORMANCE Performance::stop_region (%d,%s) %s:%d\n",CkMyPe(),
	    id_region,region_name_[id_region].c_str(),file.c_str(),line);
#endif

  int index_region = id_region;

  if (region_started_[index_region]) {

    region_started_[index_region] = false;

  } else if (warnings_) {
    if (file == "") {
      WARNING1 ("Performance::stop_region",
		"Region %s already stopped",
		region_name_[id_region].c_str());
    } else {
      WARNING3 ("Performance::stop_region",
		"Region %s already stopped %s %d",
		region_name_[id_region].c_str(),
		file.c_str(),line);
    }
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

