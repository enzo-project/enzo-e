// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      parallel_parallel.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Wed Oct 14 11:23:14 PDT 2009
/// @brief     Implementation of the Parallel class

#include "cello.h"
 
#include "parallel.hpp"

Parallel * Parallel::instance_; // (singleton design pattern)

Parallel * Parallel::instance() throw ()
{
  if (Parallel::instance_ == 0) {
#ifdef CONFIG_USE_MPI
    Parallel::instance_ = ParallelMpi::instance();
#else
    Parallel::instance_ = new Parallel;
#endif
  }
  return Parallel::instance_;
}
