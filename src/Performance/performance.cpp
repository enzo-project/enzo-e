//345678901234567890123456789012345678901234567890123456789012345678901234567890

/** 
 *********************************************************************
 *
 * @file      performance.cpp
 * @brief     
 * @author    James Bordner (jobordner@ucsd.edu)
 * @date      
 * @bug       
 * @note      
 *
 * DESCRIPTION 
 * 
 *    
 *
 * PACKAGES
 *
 *    
 * 
 * INCLUDES
 *  
 *    
 *
 * PUBLIC FUNCTIONS
 *  
 *    
 *
 * PRIVATE FUCTIONS
 *  
 *    
 *
 * $Id$
 *
 *********************************************************************
 */

#include "cello.h"

#include "error.hpp" 
#include "memory.hpp"
#include "performance.hpp"

Performance::Performance 
(
 size_t num_attributes,
 size_t num_counters,
 size_t num_groups,
 size_t num_regions)
  : counters_(),
    num_attributes_(num_attributes),
    attribute_names_(NULL),
    monotonic_attributes_(NULL),
    monotonic_attribute_values_(NULL),
    num_counters_(num_counters),
    counter_names_(NULL),
    num_groups_(num_groups),
    current_group_(0),
    group_names_(NULL),
    num_regions_(num_regions),
    region_names_(NULL),
    current_region_(0)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a Performance object with space reserved for the specified
 * number of attributes and the specified number of counters
 *
 *********************************************************************
 */
{
  Memory::begin_group("Performance");

  attribute_names_              = new std::string [ num_attributes_ ];
  counter_names_                = new std::string [ num_counters_ ];
  group_names_                  = new std::string [ num_groups_ ];
  region_names_                 = new std::string [ num_regions_ ];
  monotonic_attributes_         = new bool [num_attributes];
  for (size_t i=0; i<num_attributes; i++) {
    monotonic_attributes_[i]       = false;
  }
  monotonic_attribute_values_   = new int [num_attributes];
  for (size_t i=0; i<num_attributes; i++) {
    monotonic_attribute_values_[i] = 0;
  }

  // Create initial Counters object

  counters_.push_back(new Counters(num_attributes,num_counters));

  Memory::end_group("Performance");
}

Performance::~Performance()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Destroy a Performance object
 *
 *********************************************************************
 */
{

  Memory::begin_group("Performance");

  delete [] attribute_names_;
  delete [] counter_names_;
  delete [] group_names_;
  delete [] region_names_;
  delete [] monotonic_attributes_;
  delete [] monotonic_attribute_values_;

  for (size_t i=0; i<counters_.size(); i++) {
    delete counters_.at(i);
  }

  Memory::end_group("Performance");

}

void Performance::new_attribute(int         id_attribute, 
				 std::string attribute_name,
				 bool        is_monotonic,
				 int         max_value)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a new attribute
 *
 *********************************************************************
 */
{
}

int Performance::get_attribute(int id_attribute)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Return the value of an attribute
 *
 *********************************************************************
 */
{
  return 0;
}

void Performance::set_attribute(int id_attribute)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Assign a value to an attribute
 *
 *********************************************************************
 */
{
}

void Performance::new_group(int         id_group, 
			    std::string group_name)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a new group
 *
 *********************************************************************
 */
{
}

int Performance::get_group(int id_group)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Return the value of the currently active group
 *
 *********************************************************************
 */
{
  return 0;
}

void Performance::set_group(int id_group)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Assign a value to a group
 *
 *********************************************************************
 */
{
}

void Performance::begin_group(int group_id)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Define the start of a group
 *
 *********************************************************************
 */
{
  
  if ( current_group_ ){
    char message [ ERROR_MESSAGE_LENGTH ];
    sprintf (message, "Performance group started when one already active");
    WARNING_MESSAGE("Performance::begin_group",message);
  } else {

    current_group_ = group_id;

  }
  
}

void Performance::end_group(int group_id)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Define the end of a group
 *
 *********************************************************************
 */
{
}

void Performance::new_region(int         id_region, 
			     std::string region_name)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a new region
 *
 *********************************************************************
 */
{
}

int Performance::get_region(int id_region)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Return the value of the currently active region
 *
 *********************************************************************
 */
{
  return 0;
}

void Performance::set_region(int id_region)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Assign a value to a region
 *
 *********************************************************************
 */
{
}

void Performance::start_region(int region_id)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Define the start of a region
 *
 *********************************************************************
 */
{
}

void Performance::stop_region(int region_id)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Define the end of a region
 *
 *********************************************************************
 */
{
}


void Performance::new_counter(int         id_counter,
			      std::string counter_name)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a new user counter
 *
 *********************************************************************
 */
{
}

type_counter Performance::get_counter(int id_counter)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Return the value of a counter
 *
 *********************************************************************
 */
{
  return 0;
}

void Performance::set_counter(int          id_counter,
			      type_counter value)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Assign a value to a user counter
 *
 *********************************************************************
 */
{
}

void Performance::increment_counter(int          id_counter,
				    type_counter value)
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Increment a user counter
 *
 *********************************************************************
 */
{
}

void Performance::flush()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Flush data to disk
 *
 *********************************************************************
 */
{
}

