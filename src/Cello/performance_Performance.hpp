// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @brief    [\ref Performance] Interface for Performance class

#ifndef PERFORMANCE_PERFORMANCE_HPP
#define PERFORMANCE_PERFORMANCE_HPP

/// @def      type_counter
/// @brief    Counter value type

class Performance {

  /// @class    Performance
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Measuring and allow access to run-time
  /// parallel performance

public: // interface

  /// Initialize a Performance object
  Performance();

  /// Delete a Performance object
  ~Performance();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    WARNING("Performance::pup",
	    "skipping Performance");
    TRACE1 ("performance = %p",this);
    return;
    
    // NOTE: change this function whenever attributes change
    p | num_counters_total_;
    p | num_counters_user_;
    p | index_time_real_;
    p | papi_;
    p | counter_names_;
    p | num_regions_;
    p | region_names_;
    p | region_counters_;
    p | region_index_;
  }
#endif

  //--------------------------------------------------

  /// Begin collecting performance data
  void begin() throw();

  /// End collecting performance data
  void end() throw();

  //--------------------------------------------------

  /// Return the number of counters
  int num_counters() const throw() { return num_counters_total_; }

  ///  	Create a new user counter.
  int new_counter(std::string counter_name);

  ///  	Return the value of a counter.
  long long counter(int index_counter);

  ///  	Assign a value to a user counter.
  void assign_counter(int index_counter, long long value);

  ///  	Increment a user counter.
  void increment_counter(int index_counter, long long value);

  ///  	Increment a user counter.
  std::string counter_name (int index_counter)
  { return counter_names_[index_counter]; }

  //--------------------------------------------------

  /// Return number of regions
  int num_regions() const throw();

  /// Return the current region
  std::string region_name (int index_region) const throw();

  /// Return the index of the given region
  int region_index (std::string name) const throw();

  /// Add a new region, returning the id
  int new_region(std::string region) throw();

  /// Push a new region onto the stack
  void start_region(int index_region) throw();

  /// Push a new region onto the stack
  void stop_region(int index_region) throw();

  /// Clear the counters for the region
  void clear_region(int index_region) throw();

  /// Read the counters for the region
  void region_counters(int index_region, long long * counters) throw();

  //--------------------------------------------------

  /// Return the Papi object
  Papi * papi() { return &papi_; };

  //==================================================

private: // functions

  long long time_real_ () const
  {
    struct timeval tv;
    struct timezone tz;
    gettimeofday (&tv,&tz);
    return (long long )(1000000) * tv.tv_sec + tv.tv_usec;
  }

  //==================================================

private: // attributes

  /// Total Number of counters
  int num_counters_total_;

  /// Number of user counters
  int num_counters_user_;

  /// Index of the time counter
  int index_time_real_;

  /// PAPI counters, if available
  Papi papi_;

  /// Counter names
  std::vector<std::string> counter_names_;

  /// Counter values
  std::vector<long long> counter_values_;

  /// Number of regions in lists
  int num_regions_;

  /// list of region names
  std::vector<std::string> region_names_;

  /// list of counter values
  std::vector< std::vector<long long> > region_counters_;

  /// list of counter values
  std::vector< std::vector<long long> > region_counters_start_;


  /// mapping of region name to index
  std::map<std::string,int> region_index_;

};

#endif /* PERFORMANCE_PERFORMANCE_HPP */
