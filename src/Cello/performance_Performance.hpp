// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @brief    [\ref Performance] Interface for Performance class

#ifndef PERFORMANCE_PERFORMANCE_HPP
#define PERFORMANCE_PERFORMANCE_HPP

class Config;

/// @enum     type_counter
/// @brief    Counter value type
enum counter_type_enum {
  counter_type_unknown,
  counter_type_rel,
  counter_type_abs,
  counter_type_papi,
  counter_type_user,
  NUM_COUNTER_TYPES
};

enum index_enum {
  index_time_,
  index_bytes_,
  index_bytes_high_,
  index_bytes_highest_
};
  
/// @enum    perf_regions
/// @brief   region ID's for the Simulation performance object
enum perf_region {
  perf_simulation,
  perf_cycle,
  perf_initial,
  perf_adapt,
  perf_refresh,
  perf_compute,
  perf_output,
  perf_prepare,
  perf_last
};

class Performance {

  /// @class    Performance
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Measuring and allow access to run-time
  /// parallel performance

public: // interface

  /// Initialize a Performance object
  Performance(Config *);

  /// Delete a Performance object
  ~Performance();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    WARNING("Performance::pup",
	    "skipping Performance");
    return;
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    // NOTE: change this function whenever attributes change
    p | papi_;
    p | counter_name_;
    p | counter_type_;
    p | counter_values_;
    p | region_name_;
    p | region_counters_;
    p | region_started_;
    WARNING("Performance::pup",
	    "skipping Performance:region_index_");
    //    p | region_index_;
    WARNING("Performance::pup",
	    "skipping Performance:papi_counters_");
    //    p | papi_counters_
    p | warnings_;
  }

  /// Begin collecting performance data
  void begin() throw();

  /// End collecting performance data
  void end() throw();

  /// Return the number of counters
  int num_counters() const throw() 
  { return counter_name_.size(); }

  ///  	Create a new user counter.
  int new_counter(int counter_type, std::string counter_name);

  ///  	Return the value of a counter.
  long long counter(int index_counter) throw();

  ///  	Assign a value to a user counter.
  void assign_counter(int index_counter, long long value);

  ///  	Increment a user counter.
  void increment_counter(int index_counter, long long value);

  ///  	Return the given counter name
  std::string counter_name (int index_counter)
  { return counter_name_[index_counter]; }

  /// Return the type of the given counter index
  int counter_type (int index) const throw()
  { return counter_type_[index]; }

  ///  	Return the array of counter values
  int counter_values (const long long * values) const
  { values = &counter_values_[0];
    return counter_values_.size(); }

  /// Return number of regions
  int num_regions() const throw()
  {  return region_name_.size(); }

  /// Return the currently active region
  std::string region_name (int index_region) const throw()
  { return region_name_[index_region]; }

  /// Return the index of the given region
  int region_index (std::string name) const throw();

  /// Add a new region, returning the id
  void new_region(int index_region, std::string region) throw();

  /// Return whether performance monitoring is started for the region 
  bool is_region_active(int index_region) throw();

  /// Start counters for a code region
  void start_region(int index_region, std::string file="", int line=0,void * = 0) throw();

  /// Stop counters for a code region
  void stop_region(int index_region, std::string file="", int line=0, void * = 0) throw();

  /// Clear the counters for a code region
  void clear_region(int index_region) throw();

  /// Return counters for a code region
  void region_counters(int index_region, long long * counters) throw();

  /// Return whether the given region is active
  bool region_started(int index_region) const throw()
  { return region_started_[index_region]; }

  /// Return the associated Papi object
  Papi * papi() { return &papi_; };

private: // functions

  /// Refresh the array of current counter values
  void refresh_counters_() throw();

  /// Return the current time in usec
  long long time_real_ () const
  {
    struct timeval tv;
    struct timezone tz;
    gettimeofday (&tv,&tz);
    return (long long )(1000000) * tv.tv_sec + tv.tv_usec;
  }

  //==================================================

private: // attributes

  /// PAPI counters, if available
  Papi papi_;

  /// Counter names
  std::vector<std::string> counter_name_;

  /// Counter types
  std::vector<int> counter_type_;

  /// Counter values
  std::vector<long long> counter_values_;

  /// list of region names
  std::vector<std::string> region_name_;

  /// list of counter values
  std::vector< std::vector<long long> > region_counters_;

  /// list of counter values (should be bool but Charm++ requires int)
  std::vector< int > region_started_;

  /// mapping of region name to index
  std::map<const std::string,int> region_index_;

  /// Array for storing PAPI counter values
  long long * papi_counters_;

  /// Whether to output warning messages
  bool warnings_;

};

#endif /* PERFORMANCE_PERFORMANCE_HPP */
