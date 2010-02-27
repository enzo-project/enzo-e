/** 
 *********************************************************************
 *
 * @file      
 * @brief     
 * @author    
 * @date      
 * @ingroup
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

//======================================================================
// INCLUDES
//======================================================================

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <new>

#include "cello.h"

#include "error.hpp"
#include "memory.hpp"

/** 
 *********************************************************************
 *
 * @file      memory.cpp
 * @brief     Functions for dynamic memory management
 * @author    James Bordner
 * @date      Thu Sep  3 16:44:18 PDT 2009
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
// FUNCTIONS
//======================================================================
 
void Memory::initialize() throw ()
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
#ifdef CONFIG_USE_MEMORY
  curr_group_.push(0);
  is_active_ = true;
  group_names_[0] = strdup("");
  group_names_[1] = strdup("memory");
#endif
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
#ifdef CONFIG_USE_MEMORY
  int * buffer = (int *)(malloc(bytes + 2*sizeof(int)));


  if (buffer==0) throw ExceptionMemoryBadAllocate();

  if (is_active_) {
    buffer[0] = bytes;
    buffer[1] = curr_group_.top();

    ++ new_calls_[0] ;
    bytes_[0] += bytes;
    bytes_high_[0] = bytes_[0];

    memory_group_handle current = curr_group_.top();

    if (current != 0) {
      ++ new_calls_[current] ;
      bytes_[current] += bytes;
      bytes_high_[current] = bytes_[current];
    }

  } else {
    buffer[0] = 0;
    buffer[1] = 0;
  }
  return (void *)(buffer + 2);


#else
  return 0;
#endif
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
#ifdef CONFIG_USE_MEMORY
  int *buffer = (int *)(pointer) - 2;

  if (is_active_) {
    ++ delete_calls_[0] ;
    bytes_[0] -= buffer[0];

    memory_group_handle current = buffer[1];

    if (current != 0) {
      ++ delete_calls_[current] ;
      bytes_[current] -= buffer[0];
    }
  }

  free(buffer);
#endif
}

void Memory::new_group ( memory_group_handle group_id, const char * group_name ) throw ()
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
#ifdef CONFIG_USE_MEMORY
  if (group_id == 0 || group_id > MEMORY_MAX_NUM_GROUPS) {

    WARNING_MESSAGE("Memory::new_group()","group_id out of range");

  } else {

    group_names_[group_id] = strdup(group_name);

  }
#endif
}

void  Memory::begin_group ( memory_group_handle group_id ) throw ()
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
#ifdef CONFIG_USE_MEMORY
  // check for whether we have already called begin_group before

  bool in_range = (group_id <= MEMORY_MAX_NUM_GROUPS);

  if ( in_range ) {

    curr_group_.push(group_id);

  } else { // curr_group_ out of range

    char warning_message [ ERROR_MESSAGE_LENGTH ];

    sprintf (warning_message, "Group %d is out of range [1,%d]\n",
	     group_id, MEMORY_MAX_NUM_GROUPS);

    WARNING_MESSAGE("Memory::begin_group()",warning_message);

  }
#endif
}

void Memory::end_group ( memory_group_handle group_id ) throw ()
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
#ifdef CONFIG_USE_MEMORY

  char warning_message [ ERROR_MESSAGE_LENGTH ];

  bool in_range = (group_id <= MEMORY_MAX_NUM_GROUPS);

  if ( in_range ) {

    if (curr_group_.size() > 0) {

      if (curr_group_.top() != group_id) {
	sprintf (warning_message, 
		 "Mismatch between end_group(%d) and group stack top %d\n",
		 group_id,curr_group_.top());

	WARNING_MESSAGE("Memory::end_group",warning_message);
      }

      curr_group_.pop();

    } else {

      sprintf (warning_message, 
	       "end_group(%d) called with empty group stack\n",
	       group_id);

      WARNING_MESSAGE("Memory::end_group",warning_message);
      
    }

  } else { // curr_group_ out of range

    sprintf (warning_message, "Group %d is out of range [1,%d]\n",
	     group_id, MEMORY_MAX_NUM_GROUPS);

    WARNING_MESSAGE("Memory::end_group",warning_message);
  }
#endif
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
#ifdef CONFIG_USE_MEMORY

  // Create "null" group name if needed

  memory_group_handle current = curr_group_.top();

  if (group_names_[current] == 0) {
    group_names_[current] = strdup("");
  }

  return group_names_[current];
#else
  return "";
#endif
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
#ifdef CONFIG_USE_MEMORY
  return curr_group_.top();
#else
  return 0;
#endif
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
#ifdef CONFIG_USE_MEMORY
  return bytes_[group_handle];
#else
  return 0;
#endif
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
#ifdef CONFIG_USE_MEMORY
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::available()","");
  return 0;
#else
  return 0;
#endif
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
#ifdef CONFIG_USE_MEMORY
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::efficiency()","");
  return 0.0;
#else
  return 0;
#endif
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
#ifdef CONFIG_USE_MEMORY
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::highest(group)","");
  return 0;
#else
  return 0;
#endif
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
#ifdef CONFIG_USE_MEMORY
  check_handle_(group_handle);
  available_[group_handle] = size;
#endif
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
#ifdef CONFIG_USE_MEMORY
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::num_new()","");
  return 0;
#else
  return 0;
#endif
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
#ifdef CONFIG_USE_MEMORY
  check_handle_(group_handle);
  INCOMPLETE_MESSAGE("Memory::num_delete()","");
  return 0;
#else
  return 0;
#endif
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
#ifdef CONFIG_USE_MEMORY
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
#endif
}

void Memory::reset() throw()
{
#ifdef CONFIG_USE_MEMORY
  while (! curr_group_.empty()) curr_group_.pop();

  curr_group_.push(0);

  for (int i=0; i<MEMORY_MAX_NUM_GROUPS + 1; i++) {
    bytes_       [i] = 0;
    bytes_high_  [i] = 0;
    new_calls_   [i] = 0;
    delete_calls_[i] = 0;
  }
#endif
}

//======================================================================

#ifdef CONFIG_USE_MEMORY

bool Memory::is_active_  =  false;

std::stack<memory_group_handle> Memory::curr_group_;

char *    Memory::group_names_ [MEMORY_MAX_NUM_GROUPS + 1] = {0};

long long Memory::available_   [MEMORY_MAX_NUM_GROUPS + 1] = {0};
long long Memory::bytes_       [MEMORY_MAX_NUM_GROUPS + 1] = {0};
long long Memory::bytes_high_  [MEMORY_MAX_NUM_GROUPS + 1] = {0};
long long Memory::new_calls_   [MEMORY_MAX_NUM_GROUPS + 1] = {0};
long long Memory::delete_calls_[MEMORY_MAX_NUM_GROUPS + 1] = {0};

#endif

//======================================================================
    
#ifdef CONFIG_USE_MEMORY
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
#endif /* CONFIG_USE_MEMORY */
