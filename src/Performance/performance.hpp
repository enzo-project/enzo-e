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
  Performance(int num_attributes, 
	      int num_counters,
	      int num_groups,
	      int num_regions);

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

  /// Return the number of attributes
  size_t num_attributes();

  //--------------------------------------------------
  // GROUPS
  //--------------------------------------------------

  void new_group(int         id_group, 
		 std::string group_name);

  /// Return the value of an group
  int get_group(int id_group);

  /// Assign a value to an group
  void set_group(int id_group);

  /// Return the number of groups
  size_t num_groups();

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

  /// Return the number of regions
  size_t num_regions();

  ///  	Define the start of a region
  void start_region(int region_name);

  ///  	Define the end of a region
  void stop_region(int region_name);

  //--------------------------------------------------
  // COUNTERS
  //--------------------------------------------------

  ///  	Create a new user counter
  void new_counter(int id_counter,
		      std::string counter_name);

  ///  	Return the value of a counter
  long long get_counter(int id_counter);

  ///  	Assign a value to a user counter
  void set_counter(int id_counter,
			 long long value);
  ///  	Increment a user counter
  void increment_counter(int id_counter,
			 long long value);

  /// Return the number of counters
  size_t num_counters();

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

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Array of counters for regions
  std::vector<Counters *> counters_;

  /// Whether an attribute is monotonic
  bool * is_monotonic_;

  /// Current group
  int current_group_;

  /// Current region
  int current_region_;

  /// Whether a group is active
  bool in_group_;

  /// Whether a region is active
  bool in_region_;



};

#endif /* PERFORMANCE_HPP */
