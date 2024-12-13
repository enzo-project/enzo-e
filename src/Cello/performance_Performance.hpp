// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @brief    [\ref Performance] Interface for Performance class

#ifndef PERFORMANCE_PERFORMANCE_HPP
#define PERFORMANCE_PERFORMANCE_HPP

class Config;

class Performance {

  /// @class    Performance
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Measuring and allow access to run-time
  /// parallel performance

public: // interface

  Performance()
    :
#ifdef CONFIG_USE_PAPI
     papi_(),
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
     warnings_(false),
     index_region_current_(perf_unknown)
#ifdef CONFIG_USE_PAPI
     ,
     papi_(),
     papi_counters_(0),
#endif
#ifdef CONFIG_USE_PROJECTIONS
     ,
     projections_tracing_(),
     projections_schedule_on_(),
     projections_schedule_off_()
#endif
  { /* ... */ }

  /// Initialize a Performance object
  Performance(Config *);

  /// Delete a Performance object
  ~Performance();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    TRACEPUP;

    // NOTE: change this function whenever attributes change

#ifdef CONFIG_USE_PAPI
    p | papi_;
#endif

    p | counter_name_;
    p | counter_type_;
    p | counter_values_;
    p | counter_values_reduced_;
    p | region_name_;
    p | region_counters_;
    p | region_index_;
    p | region_multiplicity_;
    p | region_in_charm_;

#ifdef CONFIG_USE_PROJECTIONS
    p | projections_tracing_;
    p | projections_schedule_on_;
    p | projections_schedule_off_;
#endif

    p | warnings_;
    p | index_region_current_;
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

  /// Return number of regions
  int num_regions() const throw()
  {  return region_name_.size(); }

  /// Return the currently active region
  std::string region_name (int index_region) const throw()
  { return region_name_[index_region]; }

  /// Return the index of the given region
  int region_index (std::string name) const throw();

  /// Return the multiplicity (#start - #stop) of the region
  int region_multiplicity (int index_region) const throw()
  { return region_multiplicity_[index_region]; }

  /// Return whether the code region is outside the scope of Cello
  bool region_in_charm (std::string name) const throw();

  /// Add a new region, returning the id
  void new_region(int index_region, std::string region, bool in_charm=false) throw();

  /// Return whether performance monitoring is started for the region
  bool is_region_active(int index_region) throw();

  /// Start counters for a code region
  void start_region(int index_region, std::string file="", int line=0) throw();

  /// Stop counters for a code region
  void stop_region(int index_region,  std::string file="", int line=0) throw();

  /// Clear the counters for a code region
  void clear_region(int index_region) throw();

  /// Return counters for a code region
  void region_counters(int index_region, long long * counters) throw();

#ifdef CONFIG_USE_PAPI
  /// Return the associated Papi object
  Papi * papi() { return &papi_; };
#endif

#ifdef CONFIG_USE_PROJECTIONS
  /// Set whether performance tracing with projections is enabled or not
  void set_projections_tracing (bool value)
  { projections_tracing_ = value; }

  bool projections_tracing() const
  { return projections_tracing_; }

  Schedule * projections_schedule_on() const
  { return projections_schedule_on_; }

  Schedule * projections_schedule_off() const
  { return projections_schedule_off_; }
#endif

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

  /// Counter names
  std::vector<std::string> counter_name_;

  /// Counter types
  std::vector<int> counter_type_;

  /// Counter values
  std::vector<long long> counter_values_;

  /// Reduced counter values (e.g. sum over processes)
  std::vector<long long> counter_values_reduced_;

  /// list of region names
  std::vector<std::string> region_name_;

  /// list of counter values
  std::vector< std::vector<long long> > region_counters_;

  /// list of counter values (should be bool but Charm++ requires int)
  std::vector< int > region_started_;

  /// mapping of region name to index
  std::map<std::string,int> region_index_;

  /// region number of starts - number of stops
  std::vector <int> region_multiplicity_;

  /// which regions are outside scope of Cello
  std::vector<char> region_in_charm_;

  /// Whether to output warning messages
  bool warnings_;

  /// Last region index started
  int index_region_current_;

#ifdef CONFIG_USE_PAPI
  /// PAPI counters, if available
  Papi papi_;

  /// Array for storing PAPI counter values
  long long * papi_counters_;
#endif

#ifdef CONFIG_USE_PROJECTIONS
  /// Schedule for projections on / off
  bool projections_tracing_;
  Schedule * projections_schedule_on_;
  Schedule * projections_schedule_off_;
#endif
};

#endif /* PERFORMANCE_PERFORMANCE_HPP */
