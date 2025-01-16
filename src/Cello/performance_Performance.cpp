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
  counter_values_reduced_(),
  region_name_(),
  region_counters_(),
  region_index_(),
  region_multiplicity_(),
  region_in_charm_(),
#ifdef CONFIG_USE_PAPI
  papi_counters_(0),
#endif
  #ifdef CONFIG_USE_PROJECTIONS
  projections_tracing_(true),
  projections_schedule_on_(NULL),
  projections_schedule_off_(NULL),
#endif
  warnings_(config ? config->performance_warnings : false),
  index_region_current_(perf_unknown)
{

  const int in = cello::index_static();

  time_start[in] = time_real_();

  // ORDER MUST MATCH index_enum
  new_counter(counter_type_rel,"time-usec");
  // MEMORY
  new_counter(counter_type_abs,"memory_bytes-curr");
  new_counter(counter_type_abs,"memory_bytes-high");
  new_counter(counter_type_abs,"memory_bytes-highest");
  new_counter(counter_type_abs,"memory_bytes-available");

#ifdef CONFIG_USE_PAPI
  papi_.init();
#endif
#ifdef CONFIG_USE_PROJECTIONS
  int index_on = config->performance_on_schedule_index;
  int index_off = config->performance_off_schedule_index;
  projections_tracing_ = config->performance_projections_on_at_start;
  if (projections_tracing_ == false) {
    traceEnd();
  }
  if (index_on >= 0) {
    projections_schedule_on_ = Schedule::create
      ( config->schedule_var[index_on],
        config->schedule_type[index_on],
        config->schedule_start[index_on],
        config->schedule_stop[index_on],
        config->schedule_step[index_on],
        config->schedule_list[index_on]);
  }
  if (index_off >= 0) {
    projections_schedule_off_ = Schedule::create
      ( config->schedule_var[index_off],
        config->schedule_type[index_off],
        config->schedule_start[index_off],
        config->schedule_stop[index_off],
        config->schedule_step[index_off],
        config->schedule_list[index_off]);
  }
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
  counter_values_reduced_.push_back(0);

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

  counter_values_[perf_index_time]          = time_real_()-time_start[in];
  // MEMORY
  counter_values_[perf_index_bytes]         = memory->bytes();
  counter_values_[perf_index_bytes_high]    = memory->bytes_high();
  counter_values_[perf_index_bytes_highest] = memory->bytes_highest();
  counter_values_[perf_index_bytes_available] = memory->bytes_available();

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
  auto it=region_index_.find(name);
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
    region_name_.        resize(region_index+1);
    region_in_charm_.    resize(region_index+1);
    region_counters_.    resize(region_index+1);
    region_multiplicity_.resize(region_index+1);
  }

  region_index_       [region_name]    = region_index;

  region_name_        [region_index]    = region_name;
  region_in_charm_    [region_index] = in_charm;
  region_multiplicity_[region_index] = 0;

  std::vector <long long> counters;
  region_counters_    [region_index] = counters;
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

  ++region_multiplicity_[index_region];

  if (region_multiplicity_[index_region] == 1) {

    refresh_counters_();

    for (int i=0; i<num_counters(); i++) {
      if ( counter_type(i) == counter_type_abs ) {
        region_counters_[index_region][i] = counter_values_[i];
      } else {
        region_counters_[index_region][i] -= counter_values_[i];
      }
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

  --region_multiplicity_[index_region];

  if (region_multiplicity_[index_region] < 0) {
    WARNING2 ("Performance::stop_region",
              "region_multiplicity for region %s is negative %d",
              region_name(id_region).c_str(),region_multiplicity_[index_region]);
  }

  // Update region counters if final stop
  if (region_multiplicity_[index_region] == 0) {

    refresh_counters_();

    for (int i=0; i<num_counters(); i++) {
      if ( counter_type(i) == counter_type_abs ) {
        region_counters_[index_region][i] = counter_values_[i];
      } else {
        region_counters_[index_region][i] += counter_values_[i];
      }
    }
  }
}

//----------------------------------------------------------------------

bool
Performance::is_region_active(int index_region) throw()
{
  return (region_multiplicity_[index_region] > 0);
}

//----------------------------------------------------------------------

void
Performance::region_counters(int index_region, long long * counters) throw()
{
  if ( region_multiplicity_[index_region] == 0) {
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

