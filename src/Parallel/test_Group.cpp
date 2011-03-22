// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Group.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      Crashes in Parallel::initialize() in MPI_Init with LAM MPI
/// @date     Tue Apr 20 14:19:04 PDT 2010
/// @brief    Program implementing unit tests for the Group hierarchy of classes

#include "test.hpp"

#include "parallel.hpp"

int main(int argc, char ** argv)
{
  Mpi::init(int argc, char ** argv);

  Group * group_process_mpi = new GroupProcessMpi();

  unit_init(parallel->process_rank(),parallel->process_count());

  int process_count = parallel->process_count();
  int process_rank  = parallel->process_rank();
  int thread_rank   = parallel->thread_rank();

  unit_func("GroupProcessMpi","process_count");
  unit_assert(process_count == 4);
  
  if (process_count != 4) {
    parallel->finalize();
    exit(1);
  }

  unit_func("Group","Group");

  Group * group_root = new Group (0,0);
  unit_assert(group_root != NULL);
  Group * group_this = new Group (process_rank,thread_rank);
  unit_assert(group_this != NULL);

  unit_func("Group","operator ==");
  
  unit_assert ( (*grouproot == *groupthis) ==
		(process_rank == 0) );

  unit_func("Group","process_rank");
  unit_assert(unit_incomplete);

  unit_func("Group","thread_rank");
  unit_assert(unit_incomplete);

  unit_func("Group","processes");
  unit_assert(unit_incomplete);

  unit_func("Group","threads");
  unit_assert(unit_incomplete);

  unit_func("Group","group");
  unit_assert(unit_incomplete);

  unit_finalize();
  parallel->finalize();

  Mpi::finalize();
}
