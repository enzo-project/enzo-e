// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @brief    [\ref Performance] Interface for Performance class

#ifndef PERFORMANCE_PERFORMANCE_HPP
#define PERFORMANCE_PERFORMANCE_HPP

/// @enum     type_counter
/// @brief    Counter value type
#define NUM_COUNTER_TYPES 5
enum counter_type {
  counter_type_unknown,
  counter_type_basic_rel,
  counter_type_basic_abs,
  counter_type_papi,
  counter_type_user
};


enum counter_base {
  base_user      = 1000,
  base_papi      = 2000,
  base_basic_abs = 3000,
  base_basic_rel = 4000
};


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
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    WARNING("Performance::pup",
	    "skipping Performance");
    return;
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    // NOTE: change this function whenever attributes change
    p | papi_;
    PUParray (p,counter_names_,NUM_COUNTER_TYPES);
    p | counter_values_;
    p | num_regions_;
    p | region_names_;
    p | region_counters_;
    p | region_started_;
    //    p | region_index_;
    //    p | papi_counters_
    p | i0_basic_rel_;
    p | i0_basic_abs_;
    p | i0_user_;
    p | i0_papi_;
    p | n_basic_rel_;
    p | n_basic_abs_;
    p | n_user_;
    p | n_papi_;
  }
#endif

  //--------------------------------------------------

  /// Begin collecting performance data
  void begin() throw();

  /// End collecting performance data
  void end() throw();

  //--------------------------------------------------

  /// Return the number of counters
  int num_counters() const throw() { return n_basic_rel_ + n_basic_abs_ + n_papi_ + n_user_; }

  ///  	Create a new user counter.
  int new_counter(counter_type type, std::string counter_name);

  ///  	Return the value of a counter.
  long long counter(int index_counter) throw();

  ///  	Assign a value to a user counter.
  void assign_counter(int index_counter, long long value);

  ///  	Increment a user counter.
  void increment_counter(int index_counter, long long value);

  ///  	Return the given counter name
  std::string counter_name (int id)
  { counter_type type;
    int index;
    id_to_type_index_ (id,&type,&index);
    return counter_names_[type][index]; 
  }

  ///  	Return the array of counter values
  int counter_values (const long long * values) const
  { values = &counter_values_[0];
    return counter_values_.size(); }

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

  /// Return whether the given region is active
  bool region_started(int index_region) const throw()
  { return region_started_[index_region]; }

  //--------------------------------------------------

  /// Return the Papi object
  Papi * papi() { return &papi_; };

  //==================================================

  /// Convert counter ID to global index
  int id_to_index(int id_counter) const throw();

  /// Convert global index to counter ID
  int index_to_id(int index_counter) const throw();

  /// Convert counter type and type-index to counter ID
  int type_index_to_id(counter_type type, int index_type_counter) const throw();

private: // functions



  /// Get type and index relative to type from id
  void id_to_type_index_(int id_counter, counter_type * type, int * index_local) const throw();

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
  std::vector<std::string> counter_names_[NUM_COUNTER_TYPES];

  /// Counter values
  std::vector<long long> counter_values_;

  /// Number of regions in lists
  int num_regions_;

  /// list of region names
  std::vector<std::string> region_names_;

  /// list of counter values
  std::vector< std::vector<long long> > region_counters_;

  /// list of counter values (should be bool but Charm++ requires int)
  std::vector< int > region_started_;

  /// mapping of region name to index
  std::map<const std::string,int> region_index_;

  /// Array for storing PAPI counter values
  long long * papi_counters_;

  /// Indices for marking beginning of counters of different types
  int i0_basic_rel_;
  int i0_basic_abs_;
  int i0_papi_;
  int i0_user_;

  /// Indices for marking ending of counters of different types
  int n_basic_rel_;
  int n_basic_abs_;
  int n_papi_;
  int n_user_;

};

#endif /* PERFORMANCE_PERFORMANCE_HPP */
