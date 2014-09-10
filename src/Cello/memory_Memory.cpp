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

void Memory::initialize_() throw ()
{
#ifdef CONFIG_USE_MEMORY

  is_active_ = false;

  index_group_ = 0;

  if (group_name_.size() == 0) {
    new_group ("Cello");
  }

  fill_new_    = 0xaa;
  fill_delete_ = 0xdd;

  is_active_ = true;

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

    buffer[1] = index_group_;

    if (fill_new_) {
      memset (&buffer[2],fill_new_,bytes);
    }

    ++ new_calls_[0] ;
    bytes_curr_[0] += bytes;
    bytes_high_[0]    = MAX(bytes_high_[0],   bytes_curr_[0]);
    bytes_highest_[0] = MAX(bytes_highest_[0],bytes_curr_[0]);

    if (index_group_ != 0) {
      ++ new_calls_[index_group_] ;
      bytes_curr_[index_group_] += bytes;
      bytes_high_[index_group_]    = MAX(bytes_high_[index_group_],
					 bytes_curr_[index_group_]);
      bytes_highest_[index_group_] = MAX(bytes_highest_[index_group_],
					 bytes_curr_[index_group_]);
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

    int bytes = buffer[0];

    ++ delete_calls_[0] ;
    bytes_curr_[0] -= bytes;

    int index_group = buffer[1];

    if (index_group != 0) {
      ++ delete_calls_[index_group] ;
      bytes_curr_[index_group] -= bytes;
    }

    if (fill_delete_) {
      memset (&buffer[2],fill_delete_,bytes);
    }

  }

  free(buffer);

#endif
}

//----------------------------------------------------------------------

void Memory::new_group ( std::string group_name ) throw ()
/// @param  index_group    non-zero ID of the group
/// @param  group_name  Name of the group
{
#ifdef CONFIG_USE_MEMORY

  group_name_.push_back(group_name);
  bytes_limit_  .push_back(0);
  bytes_curr_   .push_back(0);
  bytes_high_   .push_back(0);
  bytes_highest_.push_back(0);
  new_calls_    .push_back(0);
  delete_calls_ .push_back(0);
  
#endif
}

//----------------------------------------------------------------------

long long Memory::bytes ( std::string group_name ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return bytes_curr_[index_group(group_name)];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

long long Memory::bytes_available ( std::string group_name ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  int index_group = this->index_group(group_name);
  if (bytes_limit_[index_group] != 0) {
    return bytes_limit_[index_group] - bytes_curr_[index_group];
  } else {
    return 0;
  }
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

float Memory::efficiency ( std::string group_name ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  int index_group = this->index_group(group_name);
  printf ("bytes_limit_[%d] = %lld\n",index_group,bytes_limit_[index_group]);
  printf ("bytes_curr_[%d] = %lld\n",index_group,bytes_curr_[index_group]);
  if (bytes_limit_[index_group] != 0) {
    return (float) bytes_curr_[index_group] / bytes_limit_[index_group];
  } else {
    return 0.0;
  }

#else
  return 0.0;
#endif
}

//----------------------------------------------------------------------

long long Memory::bytes_high ( std::string group_name ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  int index_group = this->index_group(group_name);
  TRACE1("bytes_high = %lld",bytes_high_[index_group]);
  return bytes_high_[index_group];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

long long Memory::bytes_highest ( std::string group_name ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  int index_group = this->index_group(group_name);
  TRACE1("bytes_highest = %lld",bytes_highest_[index_group]);
  return bytes_highest_[index_group];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

void Memory::set_bytes_limit ( long long size, std::string group_name )
  throw ()
{
#ifdef CONFIG_USE_MEMORY
  int index_group = this->index_group(group_name);
  bytes_limit_[index_group] = size;
#endif
}

//----------------------------------------------------------------------

long long Memory::bytes_limit ( std::string group_name ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return bytes_limit_[index_group(group_name)];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

int Memory::num_new ( std::string group_name ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return new_calls_[index_group(group_name)];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

int Memory::num_delete ( std::string group_name ) throw ()
{
#ifdef CONFIG_USE_MEMORY
  return delete_calls_[index_group(group_name)];
#else
  return 0;
#endif
}

//----------------------------------------------------------------------

void Memory::print () throw ()
{
#ifdef CONFIG_USE_MEMORY
  for (int i=0; i< group_name_.size(); i++) {
    Monitor * monitor = Monitor::instance();
    if (i == 0 || group_name_[i] != "") {
      monitor->print ("Memory","Group %s",i ? group_name_[i].c_str(): "Total");
      monitor->print ("Memory","  bytes         = %ld",long(bytes_curr_[i]));
      monitor->print ("Memory","  bytes_high    = %ld",long(bytes_high_[i]));
      monitor->print ("Memory","  bytes_highest = %ld",long(bytes_highest_[i]));
      monitor->print ("Memory","  bytes_limit   = %ld",long(bytes_limit_[i]));
      monitor->print ("Memory","  new_calls     = %ld",long(new_calls_[i]));
      monitor->print ("Memory","  delete_calls  = %ld",long(delete_calls_[i]));
    }
  }
#endif
}

//----------------------------------------------------------------------

void Memory::reset() throw()
{
#ifdef CONFIG_USE_MEMORY
  index_group_ = 0;

  for (int i=0; i<bytes_curr_.size(); i++) {
    bytes_curr_    [i] = 0;
    bytes_high_    [i] = 0;
    bytes_highest_ [i] = 0;
    new_calls_     [i] = 0;
    delete_calls_  [i] = 0;
  }
#endif
}

//----------------------------------------------------------------------

void Memory::reset_high() throw()
{
#ifdef CONFIG_USE_MEMORY
  TRACE("reset_high");
  for (int i=0; i<bytes_high_.size(); i++) {
    bytes_high_ [i] = bytes_curr_[i];
  }
#endif
}

//======================================================================

#ifdef CONFIG_USE_MEMORY

void *operator new (size_t bytes) throw (std::bad_alloc)
{
  size_t p = (size_t) Memory::instance()->allocate(bytes);
  return (void *) p;
}

#endif /* CONFIG_USE_MEMORY */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_MEMORY

void *operator new [] (size_t bytes) throw (std::bad_alloc)
{
  size_t p = (size_t) Memory::instance()->allocate(bytes);
  return (void *)(p);
}

#endif /* CONFIG_USE_MEMORY */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_MEMORY

void operator delete (void *p) throw ()
{
  if (p==0) return;
  Memory::instance()->deallocate(p);
}

#endif /* CONFIG_USE_MEMORY */

//----------------------------------------------------------------------

#ifdef CONFIG_USE_MEMORY

void operator delete [] (void *p) throw ()
{
  if (p==0) return;
  Memory::instance()->deallocate(p);
}

#endif /* CONFIG_USE_MEMORY */
