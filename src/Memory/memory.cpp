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

  if (is_active_) {
    ++ new_calls_[0] ;
    bytes_[0] += bytes;
    bytes_high_[0] = bytes_[0];

    if (curr_group_ != 0) {
      ++ new_calls_[curr_group_] ;
      bytes_[curr_group_] += bytes;
      bytes_high_[curr_group_] = bytes_[curr_group_];
    }
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

  unsigned curr_group = buffer[1];

  if (is_active_) {
    ++ delete_calls_[0] ;
    bytes_[0] -= buffer[0];

    if (curr_group != 0) {
      ++ delete_calls_[curr_group] ;
      bytes_[curr_group] -= buffer[0];
    }
  }

  free(buffer);
}

void Memory::new_group ( unsigned group_id, const char * group_name ) throw ()
/**
 *********************************************************************
 *
 * @param  group_id    non-zero ID of the group
 * @param  group_name  Name of the group
 *
 * Name a new group
 *
 *********************************************************************
 */
{
  if (group_id == 0 || group_id > MEMORY_MAX_NUM_GROUPS) {

    WARNING_MESSAGE("Memory::new_group()","group_id out of range");

  } else {

    group_names_[group_id] = strdup(group_name);

  }
}

void  Memory::begin_group ( unsigned group_id ) throw ()
/**
 *********************************************************************
 *
 * @param  group_id ID of the group
 *
 * Begin allocating memory associated with the specified group ID
 *
 * Error handling:
 *
 *  ( ) number of groups could exceed MEMORY_MAX_NUM_GROUPS_
 *        issue a warning
 *        fall back to "general" group
 *
 *********************************************************************
 */
{
  // check for whether we have already called begin_group before

  bool in_range = (group_id <= MEMORY_MAX_NUM_GROUPS);

  if ( in_range ) {

      curr_group_ = group_id;

  } else { // curr_group_ out of range

    char warning_message [ ERROR_MESSAGE_LENGTH ];

    sprintf (warning_message, "Group %d is out of range [1,%d]\n",
	     group_id, MEMORY_MAX_NUM_GROUPS);

    WARNING_MESSAGE("Memory::begin_group()",warning_message);

  }

}

void  Memory::end_group ( unsigned group_id ) throw ()
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
  bool in_range = (group_id <= MEMORY_MAX_NUM_GROUPS);

  if ( in_range ) {

    curr_group_ = 0;

  } else { // curr_group_ out of range

    char warning_message [ ERROR_MESSAGE_LENGTH ];

    sprintf (warning_message, "Group %d is out of range [1,%d]\n",
	     group_id, MEMORY_MAX_NUM_GROUPS);

    WARNING_MESSAGE("Memory::end_group()",warning_message);

  }
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


memory_group_handle Memory::current_handle () throw ()
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


long long Memory::bytes ( memory_group_handle group_handle ) throw ()
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

long long Memory::available ( memory_group_handle group_handle ) throw ()
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

long long Memory::highest ( memory_group_handle group_handle ) throw ()
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

void Memory::set_limit ( long long size, memory_group_handle group_handle )
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
  available_[group_handle] = size;
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

/// Print memory summary
void Memory::print () throw ()
/**
*********************************************************************
*
* @param  foo    
* @return        
*
* Print memory summary
*
*********************************************************************
*/
{
  for (memory_group_handle i=0; i<= MEMORY_MAX_NUM_GROUPS; i++) {

    if (i == 0 || group_names_[i] != NULL) {
      printf ("Group %s\n",i ? group_names_[i]: "Total");
      if (available_[i]) printf ("   available_    = %ld\n",long(available_[i]));
      printf ("   bytes_        = %ld\n",long(bytes_[i]));
      printf ("   bytes_high_   = %ld\n",long(bytes_high_[i]));
      printf ("   new_calls_    = %ld\n",long(new_calls_[i]));
      printf ("   delete_calls_ = %ld\n",long(delete_calls_[i]));
    }
  }
}

void Memory::reset() throw()
{
  curr_group_ = 0;

  for (int i=0; i<MEMORY_MAX_NUM_GROUPS + 1; i++) {
    bytes_       [i] = 0;
    bytes_high_  [i] = 0;
    new_calls_   [i] = 0;
    delete_calls_[i] = 0;
  }
}

//======================================================================


bool Memory::is_active_  =  true;

memory_group_handle Memory::curr_group_ = 0;

char *    Memory::group_names_ [MEMORY_MAX_NUM_GROUPS + 1] = {0};

long long Memory::available_   [MEMORY_MAX_NUM_GROUPS + 1] = {0};
long long Memory::bytes_       [MEMORY_MAX_NUM_GROUPS + 1] = {0};
long long Memory::bytes_high_  [MEMORY_MAX_NUM_GROUPS + 1] = {0};
long long Memory::new_calls_   [MEMORY_MAX_NUM_GROUPS + 1] = {0};
long long Memory::delete_calls_[MEMORY_MAX_NUM_GROUPS + 1] = {0};

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

