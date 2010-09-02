// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_ParallelGroup.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      Crashes in Parallel::initialize() in MPI_Init with LAM MPI
/// @date     Tue Apr 20 14:19:04 PDT 2010
/// @brief    Program implementing unit tests for the Group hierarchy of classes

#include "cello.hpp"
#include "parallel.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{
  Mpi::init(int argc, char ** argv);

  Group * group_process_mpi = new GroupProcessMpi();

  unit_init(parallel->process_rank(),parallel->process_count());

  int process_count = parallel->process_count();
  int process_rank  = parallel->process_rank();
  int thread_rank   = parallel->thread_rank();

  if (process_count != 4) {
    unit_assert(false);
    parallel->finalize();
    exit(1);
  }

  unit_class("Group");

  unit_func("Group");

  Group group_root (0,0);
  Group group_this (process_rank,thread_rank);

  unit_assert(true);

  unit_func("operator ==");
  
  unit_assert ( (group_root == group_this) ==
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

  Mpi::finalize();
}
