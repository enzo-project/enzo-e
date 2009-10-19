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
 int num_attributes,
 int num_counters,
 int num_groups,
 int num_regions)
  : counters_(),
    is_monotonic_(NULL),
    in_region_(false)
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

  is_monotonic_ = new bool [num_attributes];
  for (int i=0; i<num_attributes; i++) is_monotonic_[i] = false;

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

  for (size_t i=0; i<counters_.size(); i++) {
    delete [] counters_[i];
  }

  delete [] is_monotonic_;

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

size_t Performance::num_attributes()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Return the number of attributes
 *
 *********************************************************************
 */
{
  return 0;
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

size_t Performance::num_groups()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Return the number of groups
 *
 *********************************************************************
 */
{
  return 0;
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
  
  if ( in_group_ ){
    char message [ ERROR_MESSAGE_LENGTH ];
    sprintf (message, "Performance group started when one already active");
    WARNING_MESSAGE("Performance::begin_group",message);
  } else {

    in_group_ = true;
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

size_t Performance::num_regions()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Return the number of regions
 *
 *********************************************************************
 */
{
  return 0;
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

long long Performance::get_counter(int id_counter)
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

void Performance::set_counter(int       id_counter,
			      long long value)
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

void Performance::increment_counter(int       id_counter,
				    long long value)
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

size_t Performance::num_counters()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Return the number of counters
 *
 *********************************************************************
 */
{
  return 0;
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

