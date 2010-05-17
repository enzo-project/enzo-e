// $Id: test_affinity.cpp 1372 2010-04-08 05:36:42Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_affinity.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Apr 20 14:19:04 PDT 2010
/// @brief    Program implementing unit tests for the Affinity

#include "mpi.h" 
#include "parallel.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{
  Parallel * parallel = Parallel::instance();

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

  unit_class("Affinity");

  unit_func("Affinity");

  Affinity affinity_root (0,0);
  Affinity affinity_this (process_rank,thread_rank);

  unit_assert(true);

  unit_func("operator ==");
  
  unit_assert ( (affinity_root == affinity_this) ==
		(process_rank == 0) );

  parallel->finalize();
  unit_finalize();
}
