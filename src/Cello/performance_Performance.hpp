// See LICENSE_CELLO file for license and copyright information

/// @file     performance_Performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @brief    [\ref Performance] Interface for Performance class

#ifndef PERFORMANCE_PERFORMANCE_HPP
#define PERFORMANCE_PERFORMANCE_HPP

/// @def      type_counter
/// @brief    Counter value type
typedef unsigned long long type_counter;

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
    p | num_counters_;
    p | timer_;
    p | papi_;
    p | counter_names_;
    p | num_regions_;
    p | region_names_;
    p | region_counters_;
    p | region_index_;
  }
#endif

  //--------------------------------------------------

  /// Start timers and counters
  void start (int index_region = 0) throw ();

  /// Stop timers and counters
  void stop (int index_region = 0) throw ();

  double time () throw ()
  { return timer_.value(); }

  //--------------------------------------------------

  /// Return the Timer object
  Timer * timer() 
  { return &timer_; };

  /// Return the Papi object
  Papi * papi() 
  { return &papi_; };
    
  //--------------------------------------------------

  ///  	Create a new user counter.
  unsigned new_counter(std::string counter_name);

  ///  	Return the value of a counter.
  type_counter counter(unsigned id_counter);

  ///  	Assign a value to a user counter.
  void set_counter(unsigned id_counter,
		   type_counter value);
  ///  	Increment a user counter.
  void increment_counter(unsigned id_counter,
			 type_counter value);

  /// Return number of regions
  int num_regions() const throw();

  /// Return the current region
  std::string region_name (int index_region) const throw();

  /// Return the index of the given region
  int region_index (std::string name) const throw();

  /// Add a new region, returning the id
  int add_region(std::string region) throw();

  /// Push a new region onto the stack
  void start_region(int index_region) throw();

  /// Push a new region onto the stack
  void stop_region(int index_region) throw();

  /// Clear the counters for the region
  void clear_region(int index_region) throw();

  /// Read the counters for the region
  void read_region(int index_region) throw();

private: // functions

  void deallocate_ () throw ();

  type_counter time_real_ () const
  {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv,&tz);
    return (type_counter) (1000000) * tv.tv_sec + tv.tv_usec;
  }


  //--------------------------------------------------

  type_counter time_virt_ () const
  {
# ifdef CONFIG_USE_PAPI
    return PAPI_get_virt_usec();
# else
    return 0;
# endif
  }

  //--------------------------------------------------

  void new_item_ 
  (
   std::vector<std::string> & item_names,
   unsigned                 id_item, 
   std::string              item_name
   );


private: // functions

  void insert_region_(std::string) throw();

private: // attributes

  /// Number of counters
  int num_counters_;

  /// Global timer
  Timer timer_;

  /// PAPI counters, if available
  Papi papi_;

  /// Counter names
  std::vector<std::string> counter_names_;

  /// Number of regions in lists
  int num_regions_;

  /// list of region names
  std::vector<std::string> region_names_;

  /// list of counter values
  std::vector< std::vector<long long> > region_counters_;


  /// mapping of region name to index
  std::map<std::string,int> region_index_;

};

#endif /* PERFORMANCE_PERFORMANCE_HPP */
