// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file      mpi.cpp
/// @author    James Bordner (jobordner@ucsd.edu)
/// @date      Thu Oct 15 10:41:28 PDT 2009
/// @brief     MPI helper functions

#include <mpi.h>

#include "cello.h"

#include <mpi.h>
#include "parallel.hpp"

ParallelMpi ParallelMpi::instance_; // (singleton design pattern)

void ParallelMpi::initialize(int * argc, char ***argv)
{
  MPI_Init(argc,argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &size_);
  set_initialized_(true);
};

void ParallelMpi::finalize()
{
  MPI_Finalize();
  set_initialized_(false);
}
