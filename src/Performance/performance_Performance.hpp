// $Id$
// See LICENSE_CELLO file for license and copyright information

#ifndef PERFORMANCE_PERFORMANCE_HPP
#define PERFORMANCE_PERFORMANCE_HPP

/// @file     performance_Performance.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Oct 14 23:40:13 PDT 2009
/// @todo     Complete detailed description of Performance class
/// @brief    [\ref Performance] Interface for Performance class

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

enum attribute_type_enum {
  attribute_type_default,
  attribute_type_monotonic
};

class Performance {

  /// @class    Performance
  /// @ingroup  Performance
  /// @brief    [\ref Performance] Measuring and allow access to run-time
  /// parallel performance
  ///
  /// Performance data is organized into attributes, counters, groups, and
  /// regions. @@@

public: // interface

  /// Initialize a Performance object
  Performance(unsigned num_attributes, 
	      unsigned num_counters,
	      unsigned num_groups,
	      unsigned num_regions);

  /// Delete a Performance object
  ~Performance();

  /// Copy constructor
  Performance(const Performance & classname) throw();

  /// Assignment operator
  Performance & operator= (const Performance & classname) throw();

  //--------------------------------------------------

  ///  	Create a new attribute
  void new_attribute(unsigned    id_attribute, 
		     std::string attribute_name,
		     attribute_type_enum type = attribute_type_default);

  /// Return the value of an attribute
  int attribute(unsigned id_attribute);

  /// Assign a value to an attribute
  void set_attribute(unsigned id_attribute);

  //--------------------------------------------------

  void new_group(unsigned         id_group, 
		 std::string group_name);

  /// Return the value of an group
  int group(unsigned id_group);

  /// Assign a value to an group
  void set_group(unsigned id_group);

  ///  	 Define the start of a group
  void begin_group(unsigned id_group);

  ///  	Define the end of a group
  void end_group(unsigned id_group);

  //--------------------------------------------------

  void new_region(unsigned    id_region, 
		  std::string region_name);

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
  void new_counter(unsigned id_counter,
		   std::string counter_name);

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

  /// Array of counters for regions
  std::vector<Counters *> counters_;

  /// Number of attributes
  unsigned num_attributes_;

  /// Attribute names
  std::vector<std::string> attribute_names_;

  /// Which attributes are monotonic
  std::vector<attribute_type_enum> attribute_types_;

  /// Number of counters
  unsigned num_counters_;

  /// Counter names
  std::vector<std::string> counter_names_;

  /// Number of groups
  unsigned num_groups_;

  /// Current group; 0 if none
  unsigned current_group_;

  /// Group names
  std::vector<std::string> group_names_;

  /// Region names
  std::vector<std::string> region_names_;

  /// Current region; 0 if none
  unsigned current_region_;

};

#endif /* PERFORMANCE_PERFORMANCE_HPP */
