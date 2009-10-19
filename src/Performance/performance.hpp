#ifndef PERFORMANCE_HPP
#define PERFORMANCE_HPP

//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      performance.hpp
 * @brief     Classes for collecting and allowing access to performance data
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      Wed Oct 14 23:40:13 PDT 2009
 * @bug       
 * @note      
 *
 * Classes for collecting and allowing access to performance data
 *
 * $Id$
 *
 *********************************************************************
 */

#include <vector>

#include "cello.h"

#include "performance_timer.hpp"
#include "counters.hpp"

typedef unsigned long long type_counter;

class Performance {

/** 
 *********************************************************************
 *
 * @class     Performance
 * @brief     
 * @ingroup   
 *
 * 
 *
 *********************************************************************
 */

public:

  //-------------------------------------------------------------------
  // PUBLIC OPERATIONS
  //-------------------------------------------------------------------

  /// 
  Performance(size_t num_attributes, 
	      size_t num_counters,
	      size_t num_groups,
	      size_t num_regions);

  /// 
  ~Performance();

  //--------------------------------------------------
  // ATTRIBUTES
  //--------------------------------------------------

  ///  	Create a new attribute
  void new_attribute(int         id_attribute, 
		     std::string attribute_name,
		     bool        is_monotonic    = false,
		     int         max_value       = 0);

  /// Return the value of an attribute
  int get_attribute(int id_attribute);

  /// Assign a value to an attribute
  void set_attribute(int id_attribute);

  //--------------------------------------------------
  // GROUPS
  //--------------------------------------------------

  void new_group(int         id_group, 
		 std::string group_name);

  /// Return the value of an group
  int get_group(int id_group);

  /// Assign a value to an group
  void set_group(int id_group);

  ///  	 Define the start of a group
  void begin_group(int id_group);

  ///  	Define the end of a group
  void end_group(int id_group);

  //--------------------------------------------------
  // REGIONS
  //--------------------------------------------------

  void new_region(int         id_region, 
		  std::string region_name);

  /// Return the value of an region
  int get_region(int id_region);

  /// Assign a value to an region
  void set_region(int id_region);

  ///  	Define the start of a region
  void start_region(int region_name);

  ///  	Define the end of a region
  void stop_region(int region_name);

  //--------------------------------------------------
  // COUNTERS
  //--------------------------------------------------

  ///  	Create a new user counter.
  void new_counter(int id_counter,
		   std::string counter_name);

  ///  	Return the value of a counter.
  type_counter get_counter(int id_counter);

  ///  	Assign a value to a user counter.
  void set_counter(int id_counter,
		   type_counter value);
  ///  	Increment a user counter.
  void increment_counter(int id_counter,
			 type_counter value);

  //--------------------------------------------------
  // DISK
  //--------------------------------------------------

  ///  	Flush data to disk
  void flush();


private:

  //-------------------------------------------------------------------
  // PRIVATE OPERATIONS
  //-------------------------------------------------------------------

  /// 

//----------------------------------------------------------------------

  type_counter get_real_ () const
  {
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv,&tz);
    return (type_counter) (1000000) * tv.tv_sec + tv.tv_usec;
  }

//----------------------------------------------------------------------

  type_counter get_virt_ () const
  {
# ifdef CONFIG_USE_PAPI
    return PAPI_get_virt_usec();
# else
    return 0;
# endif
  }

//----------------------------------------------------------------------

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Array of counters for regions
  std::vector<Counters *> counters_;

  /// Number of attributes
  size_t num_attributes_;

  /// Attribute names
  std::string * attribute_names_;

  /// Which attributes are monotonic
  bool * monotonic_attributes_;

  /// Values of monotonic attributes
  int  * monotonic_attribute_values_;

  /// Number of counters
  size_t num_counters_;

  /// Counter names
  std::string * counter_names_;

  /// Number of groups
  size_t num_groups_;

  /// Current group; 0 if none
  int current_group_;

  /// Group names
  std::string * group_names_;

  /// Number of regions
  size_t num_regions_;

  /// Region names
  std::string * region_names_;

  /// Current region; 0 if none
  int current_region_;


};

#endif /* PERFORMANCE_HPP */
