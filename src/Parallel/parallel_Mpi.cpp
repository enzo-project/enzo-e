// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      parallel_mpi.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Oct 15 10:41:28 PDT 2009
/// @brief     MPI helper functions

#include <sstream>

#include <mpi.h>
#include "cello.h"

#include "error.hpp"
#include "parallel.hpp"

#include <boost/thread/mutex.hpp>
boost::mutex instance_mpi_mutex;

ParallelMpi * ParallelMpi::instance_mpi_ = 0; // (singleton design pattern)

void ParallelMpi::initialize(int * argc, char ***argv)
{
  MPI_Init(argc,argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);

  // Get name_ by converting rank_ to a string using std::stringstream

  std::stringstream out;
  out << rank_;
  name_ = out.str();

  Error::instance()->set_warnings_active(rank_==0);
  Error::instance()->set_incompletes_active(rank_==0);

  set_initialized_(true);
};

//----------------------------------------------------------------------

void ParallelMpi::finalize()
{
  MPI_Finalize();

  set_initialized_(false);
}

//----------------------------------------------------------------------

void ParallelMpi::abort()
{
  MPI_Abort(MPI_COMM_WORLD,1);
}

//----------------------------------------------------------------------

void ParallelMpi::halt()
{
  finalize();
  exit (0);
}

//----------------------------------------------------------------------
ParallelMpi * ParallelMpi::instance() throw ()
{ 
  // Should be thread-safe, but inefficient
  boost::mutex::scoped_lock lock(instance_mpi_mutex);
  if (ParallelMpi::instance_mpi_ == 0) {
    ParallelMpi::instance_mpi_ = new ParallelMpi;
  }
  return ParallelMpi::instance_mpi_; 
}
