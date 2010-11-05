// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_ParallelAffinity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      Crashes in Parallel::initialize() in MPI_Init with LAM MPI
/// @date     Tue Apr 20 14:19:04 PDT 2010
/// @brief    Program implementing unit tests for the ParallelAffinity

#include "cello.hpp"
#include "mpi.h" 
#include "parallel.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{
  ParallelCreate parallel_create;
  Parallel * parallel = parallel_create.create(parallel_mpi);
  parallel->initialize(&argc,&argv);

  unit_init(parallel->process_rank(),parallel->process_count());

  int process_count = parallel->process_count();
  int process_rank  = parallel->process_rank();
  int thread_rank   = parallel->thread_rank();

  if (process_count != 4) {
    unit_assert(false);
    parallel->finalize();
    exit(1);
  }

  unit_class("ParallelAffinity");

  unit_func("ParallelAffinity");

  ParallelAffinity affinity_root (0,0);
  ParallelAffinity affinity_this (process_rank,thread_rank);

  unit_assert(true);

  unit_func("operator ==");
  
  unit_assert ( (affinity_root == affinity_this) ==
		(process_rank == 0) );

  unit_func("process_rank");
  unit_assert(false);

  unit_func("thread_rank");
  unit_assert(false);

  unit_func("processes");
  unit_assert(false);

  unit_func("threads");
  unit_assert(false);

  unit_func("group");
  unit_assert(false);

  unit_finalize();
  parallel->finalize();
}
