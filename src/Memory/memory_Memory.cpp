// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      memory_Memory.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Sep  3 16:44:18 PDT 2009
/// @brief     Functions for dynamic memory management

#include "cello.hpp"

#include "memory.hpp"

#ifdef CONFIG_USE_MEMORY
Memory Memory::instance_; // (singleton design pattern)
#endif

//======================================================================
// FUNCTIONS
//======================================================================
 
void Memory::initialize_() throw ()
{
#ifdef CONFIG_USE_MEMORY
  curr_group_.push(0);
  is_active_ = true;

  for (int i=0; i<MEMORY_MAX_NUM_GROUPS + 1; i++) {
    group_names_ [i] = 0;
    limit_       [i] = 0;
    bytes_       [i] = 0;
    bytes_high_  [i] = 0;
    new_calls_   [i] = 0;
    delete_calls_[i] = 0;
  }
  group_names_[0] = strdup("");       // MEMORY LEAK
  group_names_[1] = strdup("memory"); // MEMORY LEAK

#endif
}

//----------------------------------------------------------------------

void Memory::finalize_() throw ()
{
#ifdef CONFIG_USE_MEMORY
  for (int i=0; i<MEMORY_MAX_NUM_GROUPS + 1; i++) {
    delete [] group_names_[i];
  }
#endif
}

//----------------------------------------------------------------------

void * Memory::allocate ( size_t bytes ) throw (ExceptionMemoryBadAllocate())
/// @param  bytes   Number of bytes to allocate
/// @return        Pointer to the allocated memory
{
#ifdef CONFIG_USE_MEMORY
  int * buffer = (int *)(malloc(bytes + 2*sizeof(int)));


  if (buffer==0) throw ExceptionMemoryBadAllocate();

  if (is_active_) {
    buffer[0] = bytes;
    buffer[1] = curr_group_.top();

    ++ new_calls_[0] ;
    bytes_[0] += bytes;
    bytes_high_[0] = MAX(bytes_high_[0],bytes_[0]);

    memory_group_handle current = curr_group_.top();

    if (current != 0) {
      ++ new_calls_[current] ;
      bytes_[current] += bytes;
      bytes_high_[current] = MAX(bytes_high_[current],bytes_[current]);
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

//----------------------------------------------------------------------

void Memory::deallocate ( void * pointer )
  throw (ExceptionMemoryBadDeallocate())
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

//----------------------------------------------------------------------

void Memory::new_group ( memory_group_handle group_id, const char * group_name ) throw ()
/// @param  group_id    non-zero ID of the group
/// @param  group_name  Name of the group
{
#ifdef CONFIG_USE_MEMORY
  if (group_id == 0 || group_id > MEMORY_MAX_NUM_GROUPS) {

    WARNING("Memory::new_group()","group_id out of range");

  } else {

    group_names_[group_id] = strdup(group_name);

  }
#endif
}

//----------------------------------------------------------------------

void  Memory::begin_group ( memory_group_handle group_id ) throw ()
/// @param  group_id ID of the group
///
/// If number of groups exceeds MEMORY_MAX_NUM_GROUPS_, issue a
/// warning and fall back to "general" group
{
#ifdef CONFIG_USE_MEMORY
  // check for whether we have already called begin_group before

  bool in_range = (group_id <= MEMORY_MAX_NUM_GROUPS);

  if ( in_range ) {

    curr_group_.push(group_id);

  } else { // curr_group_ out of range

    char warning_message [ ERROR_LENGTH ];

    sprintf (warning_message, "Group %d is out of range [1,%d]\n",
	     group_id, MEMORY_MAX_NUM_GROUPS);

    WARNING("Memory::begin_group()",warning_message);

  }
#endif
}

//----------------------------------------------------------------------

void Memory::end_group ( memory_group_handle group_id ) throw ()
{
#ifdef CONFIG_USE_MEMORY

  char warning_message [ ERROR_LENGTH ];

  bool in_range = (group_id <= MEMORY_MAX_NUM_GROUPS);

  if ( in_range ) {

    if (curr_group_.size() > 0) {

      if (curr_group_.top() != group_id) {
	sprintf (warning_message, 
		 "Mismatch between end_group(%d) and group stack top %d\n",
		 group_id,curr_group_.top());

	WARNING("Memory::end_group",warning_message);
      }

      curr_group_.pop();

    } else {

      sprintf (warning_message, 
	       "end_group(%d) called with empty group stack\n",
	       group_id);

      WARNING("Memory::end_group",warning_message);
      
    }

  } else { // curr_group_ out of range

    sprintf (warning_message, "Group %d is out of range [1,%d]\n",
	     group_id, MEMORY_MAX_NUM_GROUPS);

    WARNING("Memory::end_group",warning_message);
  }
#endif
}

//----------------------------------------------------------------------

memory_group_handle Memory::current_group (const char ** group_name) throw ()
{
#ifdef CONFIG_USE_MEMORY

  // Create "null" group name if needed

  memory_group_handle current = curr_group_.top();

  if (group_names_[current] == 0) {
    group_names_[current] = strdup("");
  }

  *group_name = group_names_[current];
  return current;

#else

  *group_name = 0;
  return 0;

#endif
}

//----------------------------------------------------------------------

long long Memory::bytes ( memory_group_handle group_handle ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return bytes_[group_handle];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

long long Memory::available ( memory_group_handle group_handle ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  if (limit_[group_handle] != 0) {
    return limit_[group_handle] - bytes_[group_handle];
  } else {
    return 0;
  }
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

float Memory::efficiency ( memory_group_handle group_handle ) throw ()
{
#ifdef CONFIG_USE_MEMORY

  if (limit_[group_handle] != 0) {
    return (float) bytes_[group_handle] / limit_[group_handle];
  } else {
    return 0.0;
  }

#else
  return 0.0;
#endif
}

//----------------------------------------------------------------------

long long Memory::highest ( memory_group_handle group_handle ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return bytes_high_[group_handle];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

void Memory::set_limit ( long long size, memory_group_handle group_handle )
  throw ()
{
#ifdef CONFIG_USE_MEMORY
  limit_[group_handle] = size;
#endif
}

//----------------------------------------------------------------------

long long Memory::limit ( memory_group_handle group_handle ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return limit_[group_handle];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

int Memory::num_new ( memory_group_handle group_handle ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return new_calls_[group_handle];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

int Memory::num_delete ( memory_group_handle group_handle ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return delete_calls_[group_handle];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

void Memory::print () throw ()
{
#ifdef CONFIG_USE_MEMORY
  for (memory_group_handle i=0; i<= MEMORY_MAX_NUM_GROUPS; i++) {
    if (i == 0 || group_names_[i] != NULL) {
      PARALLEL_PRINTF ("Group %s\n",i ? group_names_[i]: "Total");
      PARALLEL_PRINTF ("   limit        = %ld\n",long(limit_[i]));
      PARALLEL_PRINTF ("   bytes        = %ld\n",long(bytes_[i]));
      PARALLEL_PRINTF ("   bytes_high   = %ld\n",long(bytes_high_[i]));
      PARALLEL_PRINTF ("   new_calls    = %ld\n",long(new_calls_[i]));
      PARALLEL_PRINTF ("   delete_calls = %ld\n",long(delete_calls_[i]));
    }
  }
#endif
}

//----------------------------------------------------------------------

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
