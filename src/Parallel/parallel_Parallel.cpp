// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      parallel_Parallel.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Oct 14 11:23:14 PDT 2009
/// @brief     Implementation of the Parallel class

#include "cello.h"
 
#include "parallel.hpp"

#include <boost/thread/mutex.hpp>
boost::mutex instance_mutex;

Parallel * Parallel::instance_ = 0; // (singleton design pattern)

Parallel * Parallel::instance() throw ()
{
  // Should be thread-safe, but inefficient
  boost::mutex::scoped_lock lock(instance_mutex);
  if (Parallel::instance_ == 0) {
#ifdef CONFIG_USE_MPI
    Parallel::instance_ = ParallelMpi::instance();
#else
    Parallel::instance_ = new Parallel;
#endif
  }
  return Parallel::instance_;
}
