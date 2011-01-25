// $Id: memory_override.cpp 1688 2010-08-03 22:34:22Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     memory_new_delete.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Nov 10 12:01:51 PST 2010
/// @brief    Override C++ new[] and delete[]

//======================================================================

#include "cello.hpp"

#include "memory.hpp"

#ifdef CONFIG_USE_MEMORY
void *operator new (size_t bytes) throw (std::bad_alloc)
{
  size_t p = (size_t) Memory::instance()->allocate(bytes);

  // Return pointer to new storage

  return (void *) p;
}

//----------------------------------------------------------------------

void *operator new [] (size_t bytes) throw (std::bad_alloc)
{
  size_t p = (size_t) Memory::instance()->allocate(bytes);

  // Return pointer to new storage

  return (void *)(p);

}

//----------------------------------------------------------------------

void operator delete (void *p) throw ()
{
  if (p==0) return;

  Memory::instance()->deallocate(p);

}

//----------------------------------------------------------------------

void operator delete [] (void *p) throw ()
{
  if (p==0) return;

  Memory::instance()->deallocate(p);
}
#endif /* CONFIG_USE_MEMORY */
