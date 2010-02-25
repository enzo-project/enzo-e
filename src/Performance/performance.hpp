//
// $Id$
//
// See LICENSE_CELLO file for license and copyright information
//

#ifndef PERFORMANCE_HPP
#define PERFORMANCE_HPP

/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
 * @bug       
 * @note      
 *
 *--------------------------------------------------------------------
 *
 * DESCRIPTION:
 *
 *    
 *
 * CLASSES:
 *
 *    
 *
 * FUCTIONS:
 *
 *    
 *
 * USAGE:
 *
 *    
 *
 * REVISION HISTORY:
 *
 *    
 *
 * COPYRIGHT: See the LICENSE_CELLO file in the project directory
 *
 *--------------------------------------------------------------------
 *
 * $Id$
 *
 *********************************************************************
 */


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

enum enum_item_type {
  item_unknown,
  item_attribute,
  item_counter,
  item_group,
  item_region
};

/** 
 *********************************************************************
 *
 * @class     Performance
 * @brief     
 * @ingroup   Performance
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
  Performance(unsigned num_attributes, 
	      unsigned num_counters,
	      unsigned num_groups,
	      unsigned num_regions);

  /// 
  ~Performance();

  //--------------------------------------------------
  // ATTRIBUTES
  //--------------------------------------------------

  ///  	Create a new attribute
  void new_attribute(unsigned    id_attribute, 
		     std::string attribute_name,
		     bool        is_monotonic    = false);

  /// Return the value of an attribute
  int get_attribute(unsigned id_attribute);

  /// Assign a value to an attribute
  void set_attribute(unsigned id_attribute);

  //--------------------------------------------------
  // GROUPS
  //--------------------------------------------------

  void new_group(unsigned         id_group, 
		 std::string group_name);

  /// Return the value of an group
  int get_group(unsigned id_group);

  /// Assign a value to an group
  void set_group(unsigned id_group);

  ///  	 Define the start of a group
  void begin_group(unsigned id_group);

  ///  	Define the end of a group
  void end_group(unsigned id_group);

  //--------------------------------------------------
  // REGIONS
  //--------------------------------------------------

  void new_region(unsigned    id_region, 
		  std::string region_name);

  /// Return the value of an region
  int get_region(unsigned id_region);

  /// Assign a value to an region
  void set_region(unsigned id_region);

  ///  	Define the start of a region
  void start_region(unsigned region_name);

  ///  	Define the end of a region
  void stop_region(unsigned region_name);

  //--------------------------------------------------
  // COUNTERS
  //--------------------------------------------------

  ///  	Create a new user counter.
  void new_counter(unsigned id_counter,
		   std::string counter_name);

  ///  	Return the value of a counter.
  type_counter get_counter(unsigned id_counter);

  ///  	Assign a value to a user counter.
  void set_counter(unsigned id_counter,
		   type_counter value);
  ///  	Increment a user counter.
  void increment_counter(unsigned id_counter,
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

  void new_item_ 
  (
   enum_item_type  item_type,
   std::string item_type_name,
   unsigned id_item, 
   std::string item_name,
   std::string * item_names,
   unsigned num_items_
   );

private:

  //-------------------------------------------------------------------
  // PRIVATE ATTRIBUTES
  //-------------------------------------------------------------------

  /// Array of counters for regions
  std::vector<Counters *> counters_;

  /// Number of attributes
  unsigned num_attributes_;

  /// Attribute names
  std::string * attribute_names_;

  /// Which attributes are monotonic
  bool * monotonic_attributes_;

  /// Values of monotonic attributes; 0 for non-monotonic
  int  * monotonic_attribute_values_;

  /// Number of counters
  unsigned num_counters_;

  /// Counter names
  std::string * counter_names_;

  /// Number of groups
  unsigned num_groups_;

  /// Current group; 0 if none
  unsigned current_group_;

  /// Group names
  std::string * group_names_;

  /// Number of regions
  unsigned num_regions_;

  /// Region names
  std::string * region_names_;

  /// Current region; 0 if none
  unsigned current_region_;


};

#endif /* PERFORMANCE_HPP */
