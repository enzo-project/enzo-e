//345678901234567890123456789012345678901234567890123456789012345678901234567890

/*
 * ENZO: THE NEXT GENERATION
 *
 * A parallel astrophysics and cosmology application
 *
 * Copyright (C) 2009 James Bordner
 * Copyright (C) 2009 Laboratory for Computational Astrophysics
 * Copyright (C) 2009 Regents of the University of California
 *
 * See CELLO_LICENSE in the main directory for full license agreement
 *
 */


/** 
 *********************************************************************
 *
 * @file      memory.cpp
 * @brief     Functions for dynamic memory management
 * @author    James Bordner
 * @date      Thu Sep  3 16:44:18 PDT 2009
 * @bug       
 * @note      
 *
 * DESCRIPTION 
 * 
 *    Functions for dynamic memory management
 *
 * PACKAGES
 *
 *    NONE
 * 
 * INCLUDES
 *  
 *    NONE
 *
 * PUBLIC FUNCTIONS
 *  
 *   ( ) static Memory();
 *   ( ) static void * allocate(size_t size);
 *   ( ) static void * allocate(size_t size, std::string class);
 *   ( ) static deallocate();
 *   ( ) static deallocate(std::string class);
 *   ( ) static current(class);
 *   ( ) static available();
 *   ( ) static efficiency();
 *   ( ) static highest(class);
 *   ( ) static set_highest();
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

//======================================================================
// INCLUDES
//======================================================================

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <new>

#include "error.hpp"
#include "memory.hpp"

//======================================================================
// FUNCTIONS
//======================================================================
 
Memory::Memory() throw ()
/**
 *********************************************************************
 *
 * @param         There are no parameters
 * @return        There is no return value
 *
 * Initialize the Memory class
 *
 *********************************************************************
 */
{
  printf ("Memory::Memory()\n");
  group_names_[0] = strdup("");
  group_names_[1] = strdup("memory");
}

void * Memory::allocate ( size_t bytes ) throw (ExceptionMemoryBadAllocate())
/**
 *********************************************************************
 *
 * @param  bytes   Number of bytes to allocate
 * @return        Pointer to the allocated memory
 *
 * Allocate memory with the default group
 *
 *********************************************************************
 */
{
  int * buffer = (int *)(malloc(bytes + 2*sizeof(int)));
  if (buffer==0) throw ExceptionMemoryBadAllocate();

  buffer[0] = bytes;
  buffer[1] = curr_group_;

  ++ new_calls_[curr_group_] ;
  new_bytes_[curr_group_] += bytes;
  bytes_[curr_group_] += bytes;
  if (bytes_[curr_group_] > bytes_high_[curr_group_]) {
    bytes_high_[curr_group_] = bytes_[curr_group_];
  }

  return (void *)(buffer + 2);
}

void Memory::deallocate ( void * pointer )
  throw (ExceptionMemoryBadDeallocate())
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  int *buffer = (int *)(pointer) - 2;

  int curr_group = buffer[1];

  ++ delete_calls_[curr_group] ;
  delete_bytes_[curr_group] += buffer[0];
  bytes_[curr_group]        -= buffer[0];

  free(buffer);
}

void  Memory::begin_group ( const char * group_name ) throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Begin allocating memory associated with the specified group name
 *
 * Error handling:
 *
 *  ( ) number of groups could exceed MEMORY_MAX_NUM_GROUPS_
 *        issue a warning
 *        fall back to "general" group
 *
 *  ( ) begin_group/end_group should be properly nested
 *        issue a warning
 *        fall back to "general" group
 *
 *********************************************************************
 */
{
  // check for whether we have already called begin_group before
  if (curr_group_ != 0) {
    INCOMPLETE_MESSAGE("Memory::begin_group()","");
    curr_group_ = 0;
  }

  // Test if group already exists, and set curr_group_ if it does

  for (int i=1; i<num_groups_; i++) {

    if (strcmp(group_names_[i],group_name) == 0) {

      curr_group_ = (memory_group_handle) (i);

    }

  }

  if (curr_group_ == 0 && (num_groups_ < MEMORY_MAX_NUM_GROUPS - 1)) {

    // Otherwise, add a new group if there is room, and set curr_group_

    group_names_[num_groups_] = strdup(group_name);
    curr_group_ = num_groups_;

    num_groups_++;

  } else {
    
    // Otherwise we hit the array limit: set curr_group_ to 0
  
    INCOMPLETE_MESSAGE("Memory::begin_group()","");
    // WARNING MESSAGE: reverting to "null" group 0

    curr_group_ = 0;
  }

}

void  Memory::end_group ( const char * group_name ) throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 *
 *********************************************************************
 */
{
  // check for whether we have already called begin_group before
  if (curr_group_ == 0) {
    INCOMPLETE_MESSAGE("Memory::end_group()","");
  }
  curr_group_ = 0;
}


const char * Memory::current_group () throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 *
 *********************************************************************
 */
{
  // Create "null" group name if needed

  if (group_names_[0] == 0) {
    group_names_[0] = strdup("");
  }

  return group_names_[curr_group_];
}


int Memory::current_handle () throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 *
 *********************************************************************
 */
{
  return curr_group_;
}


size_t Memory::bytes ( memory_group_handle group_handle ) throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  return bytes_[group_handle];
}

size_t Memory::available ( memory_group_handle group_handle ) throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::available()","");
  return 0;
}

float Memory::efficiency ( memory_group_handle group_handle ) throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::efficiency()","");
  return 0.0;
}

size_t Memory::highest ( memory_group_handle group_handle ) throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::highest(group)","");
  return 0;
}

void Memory::set_limit ( size_t size, memory_group_handle group_handle )
  throw ()
/**
 *********************************************************************
 *
 * @param  foo    Description of argument foo
 * @return        There is no return value
 *
 * Detailed description of the function
 *
 *********************************************************************
 */
{
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::set_highest()","");
  return;
}

/// Return the number of calls to allocate for the group
int Memory::num_new ( memory_group_handle group_handle ) throw ()
/**
*********************************************************************
*
* @param  foo    Description of argument foo
* @return        There is no return value
*
* Detailed description of the function
*
*********************************************************************
*/
{
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::num_new()","");
  return 0;
}

/// Return the number of calls to deallocate for the group
int Memory::num_delete ( memory_group_handle group_handle ) throw ()
/**
*********************************************************************
*
* @param  foo    Description of argument foo
* @return        There is no return value
*
* Detailed description of the function
*
*********************************************************************
*/
{
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::num_delete()","");
  return 0;
}

//======================================================================


int    Memory::num_groups_ = 1;

int    Memory::curr_group_ = 0;

char * Memory::group_names_ [MEMORY_MAX_NUM_GROUPS] = {0};

size_t Memory::bytes_[MEMORY_MAX_NUM_GROUPS] = {0};

size_t Memory::bytes_high_[MEMORY_MAX_NUM_GROUPS] = {0};

size_t Memory::new_calls_[MEMORY_MAX_NUM_GROUPS] = {0};

size_t Memory::new_bytes_[MEMORY_MAX_NUM_GROUPS] = {0};

size_t Memory::delete_calls_[MEMORY_MAX_NUM_GROUPS] = {0};

size_t Memory::delete_bytes_[MEMORY_MAX_NUM_GROUPS] = {0};

//======================================================================
    
void *operator new (size_t bytes) throw (std::bad_alloc)
{

  size_t p = (size_t) Memory::allocate(bytes);


  // Return pointer to new storage

  return (void *) p;
}

//----------------------------------------------------------------------

void *operator new [] (size_t bytes) throw (std::bad_alloc)
{

  size_t p = (size_t) Memory::allocate(bytes);


  // Return pointer to new storage

  return (void *)(p);

}

//----------------------------------------------------------------------

void operator delete (void *p) throw ()
{
  if (p==0) return;

  Memory::deallocate(p);

}

//----------------------------------------------------------------------

void operator delete [] (void *p) throw ()
{
  if (p==0) return;


  Memory::deallocate(p);

}

