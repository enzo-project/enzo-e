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
 
#include "performance.hpp"

Performance::Performance()
/**
 *********************************************************************
 *
 * @param         
 * @return        
 *
 * Create a Performance object
 *
 *********************************************************************
 */
{
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
}

void Performance::group_begin(std::string)
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
}

void Performance::group_end(std::string group_name)
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

void Performance::region_start(std::string region_name)
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

void Performance::region_stop(std::string region_name)
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


void Performance::attribute_create(type_attribute id_attribute, 
				   std::string    attribute_name,
				   bool           is_monotonic,
				   int            max_value)
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

int Performance::attribute_get(type_attribute id_attribute)
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

void Performance::attribute_set(type_attribute id_attribute)
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

size_t Performance::attribute_count()
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

void Performance::counter_create(type_counter id_counter,
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

long long Performance::counter_get(type_counter id_counter)
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

void Performance::counter_set(type_counter id_counter,
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

void Performance::counter_increment(type_counter id_counter,
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

size_t Performance::counter_count()
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

