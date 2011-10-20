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
 
Memory::~Memory() throw ()
{
  // for (int i=0; i<max_group_id_ +1; i++) {
  //   free (group_names_[i]);
  //   group_names_[i] = 0;
  // }
}

//======================================================================

void Memory::initialize_() throw ()
{
#ifdef CONFIG_USE_MEMORY
  curr_group_.push(0);
  max_group_id_ = 0;
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
  // WARNING: not called by end of program

#ifdef CONFIG_USE_MEMORY
  for (int i=0; i<MEMORY_MAX_NUM_GROUPS + 1; i++) {
    delete [] group_names_[i];
  }
#endif
}

//----------------------------------------------------------------------

void * Memory::allocate ( size_t bytes ) throw ()
/// @param  bytes   Number of bytes to allocate
/// @return        Pointer to the allocated memory
{
#ifdef CONFIG_USE_MEMORY
  int * buffer = (int *)(malloc(bytes + 2*sizeof(int)));


  ASSERT("Memory::allocate",
	 "Cannot allocate buffer: out of memory",
	 buffer);

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

void Memory::deallocate ( void * pointer ) throw()
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

  ASSERT ("Memory::new_group()","group_id out of range",
	  1 <= group_id && group_id <= MEMORY_MAX_NUM_GROUPS);

  group_names_[group_id] = strdup(group_name);

  max_group_id_  = MAX(max_group_id_, group_id);

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

  bool in_range = (group_id <= max_group_id_);

  if ( in_range ) {

    curr_group_.push(group_id);

  } else { // curr_group_ out of range

    WARNING2("Memory::begin_group()",
	     "Group %d is out of range [1,%d]\n",
	     group_id, max_group_id_);

  }
#endif
}

//----------------------------------------------------------------------

void Memory::end_group ( memory_group_handle group_id ) throw ()
{
#ifdef CONFIG_USE_MEMORY

  bool in_range = (group_id <= max_group_id_);

  if ( in_range ) {

    if (curr_group_.size() > 0) {

      if (curr_group_.top() != group_id) {

	WARNING2("Memory::end_group",
		"Mismatch between end_group(%d) and group stack top %d\n",
		 group_id,curr_group_.top());
      }

      curr_group_.pop();

    } else {

      WARNING1("Memory::end_group",
	       "end_group(%d) called with empty group stack\n",
	       group_id);
      
    }

  } else { // curr_group_ out of range

    WARNING2("Memory::end_group",
	     "Group %d is out of range [1,%d]\n",
	     group_id, max_group_id_);
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

long long Memory::bytes_high ( memory_group_handle group_handle ) throw ()
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
  for (memory_group_handle i=0; i<= max_group_id_; i++) {
    Monitor * monitor = Monitor::instance();
    if (i == 0 || group_names_[i] != NULL) {
      monitor->print ("Memory","Group %s",i ? group_names_[i]: "Total");
      monitor->print ("Memory","  limit        = %ld",long(limit_[i]));
      monitor->print ("Memory","  bytes        = %ld",long(bytes_[i]));
      monitor->print ("Memory","  bytes_high   = %ld",long(bytes_high_[i]));
      monitor->print ("Memory","  new_calls    = %ld",long(new_calls_[i]));
      monitor->print ("Memory","  delete_calls = %ld",long(delete_calls_[i]));
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

  for (int i=0; i<max_group_id_ + 1; i++) {
    bytes_       [i] = 0;
    bytes_high_  [i] = 0;
    new_calls_   [i] = 0;
    delete_calls_[i] = 0;
  }
#endif
}

//----------------------------------------------------------------------

void Memory::reset_high() throw()
{
#ifdef CONFIG_USE_MEMORY
  for (int i=0; i<max_group_id_ + 1; i++) {
    bytes_high_ [i] = bytes_[i];
  }
#endif
}

//----------------------------------------------------------------------

void Memory::check_handle_(memory_group_handle group_handle) throw ()
{  
#ifdef CONFIG_USE_MEMORY
  ASSERT2 ("Memory::check_handle_",
	   "group_handle %d is out of range [0:%d]",
	   int(group_handle), max_group_id_,
	   0 <= group_handle && group_handle <= max_group_id_);
#endif
}
