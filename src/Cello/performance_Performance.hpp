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
    p | timer_;
    p | papi_;
    WARNING("Performance::pup",
	    "skipping counters_ [ not accessed except deallocate ]");
  //  std::vector<PerfCounters *> counters_;
    p | counter_names_;
    p | region_names_;
    p | current_region_;
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

  unsigned new_region(std::string region_name);

  /// Return the value of an region
  int region(unsigned id_region);

  /// Assign a value to an region
  void set_region(unsigned id_region);

  ///  	Define the start of a region
  void start_region(unsigned region_name);

  ///  	Define the end of a region
  void stop_region(unsigned region_name);

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

  //--------------------------------------------------

  ///  	Flush data to disk
  void flush();

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

private: // attributes

  /// Global timer
  Timer timer_;

  /// PAPI counters, if available
  Papi papi_;

  /// Array of counters for regions
  std::vector<PerfCounters *> counters_;

  /// Counter names
  std::vector<std::string> counter_names_;

  /// Region names
  std::vector<std::string> region_names_;

  /// Current region; 0 if none
  unsigned current_region_;


};

#endif /* PERFORMANCE_PERFORMANCE_HPP */
