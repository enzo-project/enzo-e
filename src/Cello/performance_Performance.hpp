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

enum item_enum {
  item_unknown,
  item_attribute,
  item_counter,
  item_group,
  item_region
};

class Monitor;
class Performance {

  /// @class    Performance
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Measuring and allow access to run-time
  /// parallel performance
  ///
  /// Performance data is organized into attributes, counters, groups, and
  /// regions

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
    // NOTE: change this function whenever attributes change
    p | timer_;
    p | papi_;
    WARNING("Performance::pup",
	    "skipping std::vector<PerfCounters *> counters_");
  //  std::vector<PerfCounters *> counters_;
    p | attribute_names_;
    p | counter_names_;
    p | group_names_;
    p | region_names_;
    p | attribute_monotonic_;
    p | current_group_;
    p | current_region_;
  }
#endif

  //--------------------------------------------------

  /// Start timers and counters
  void start () throw ();

  /// Stop timers and counters
  void stop () throw ();

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

  ///  	Create a new attribute
  unsigned new_attribute( std::string attribute_name,
			  bool is_monotonic = false);

  /// Return the value of an attribute
  int attribute(unsigned id_attribute);

  /// Assign a value to an attribute
  void set_attribute(unsigned id_attribute, int value);

  //--------------------------------------------------

  unsigned new_group(std::string group_name);

  /// Return the value of an group
  int group(unsigned id_group);

  /// Assign a value to an group
  void group_set(unsigned id_group);

  ///  	 Define the start of a group
  void begin_group(unsigned id_group);

  ///  	Define the end of a group
  void end_group(unsigned id_group);

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

  // void print_rusage_ (const Monitor * monitor) const throw ();

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

  /// Attribute names
  std::vector<std::string> attribute_names_;

  /// Counter names
  std::vector<std::string> counter_names_;

  /// Group names
  std::vector<std::string> group_names_;

  /// Region names
  std::vector<std::string> region_names_;

  /// Which attributes are monotonic
  std::vector<int> attribute_monotonic_;

  /// Current group; 0 if none
  unsigned current_group_;

  /// Current region; 0 if none
  unsigned current_region_;


};

#endif /* PERFORMANCE_PERFORMANCE_HPP */
