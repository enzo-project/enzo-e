// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_GroupProcess.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      Crashes in Parallel::initialize() in MPI_Init with LAM MPI
/// @todo     Remove need for dynamic casts for MPI-specific code
/// @date     Tue Apr 20 14:19:04 PDT 2010
/// @brief    Program implementing unit tests for the Parallel class

#include "test.hpp"

#include "parallel.hpp"

//----------------------------------------------------------------------

void init_array(double * array, int length, int rank)
{
  for (int i=0; i<length-1; i++) {
    array[i] = 1.0*rank;
  }
  array[length-1] = -1.0*rank;
}

//----------------------------------------------------------------------

bool test_array(double * array, int length, int rank, int value)
{
  bool valid = true;
  for (int i=0; i<length-1; i++) {
    valid = valid && (array[i] == 1.0*value);
  }
  valid = valid && (array[length-1] == -1.0*rank);
  if (!valid) {
    char message[ERROR_LENGTH];
    sprintf (message,
	     "%d expected [%g %g %g %g %g]; actual [%g %g %g %g %g]\n",
	     rank,
	     1.0*value,1.0*value,1.0*value,1.0*value,-1.0*rank,
	     array[0], array[1], array[2], array[3], array[4]);
    WARNING("test_array",message);
  }
  return valid;
}

//======================================================================

#include PARALLEL_CHARM_INCLUDE(test_GroupProcess.decl.h)

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  int np;
  if (PARALLEL_ARGC != 2) {
    PARALLEL_PRINTF ("Usage: %s <num-procs>\n",PARALLEL_ARGV[0]);
    PARALLEL_EXIT;
  } else {
    np = atoi(PARALLEL_ARGV[1]);
  }


  GroupProcess * group_process = GroupProcess::create();

  int rank = group_process->rank();
  int size = group_process->size();

  unit_init(group_process->rank(), group_process->size());

  unit_func("GroupProcess","size");
  unit_assert(size == np);

  unit_func("GroupProcess","rank");
  unit_assert(0 <= rank && rank < np);

  unit_func("GroupProcess","is_root");
  unit_assert(group_process->is_root() == (rank == 0));

  // Test that init_array() and test_array() work independently of Parallel
  const int n = np;
  double array_source[n+1], array_dest[n+1];

  unit_func("GroupProcess","wait");

#ifdef CONFIG_USE_MPI
  const int num_blocking = 2;
#else
  const int num_blocking = 1;
#endif

  for (int blocking_send = 0; blocking_send < num_blocking; blocking_send++) {
    for (int blocking_recv = 0; blocking_recv < num_blocking; blocking_recv++) {

      init_array(array_source,n+1,rank);
      init_array(array_dest,  n+1,rank);
      unit_assert(test_array(array_source,n+1,rank,rank));

#ifdef CONFIG_USE_MPI

      // Blocking send or receive is specific to MPI

      GroupProcessMpi * group_process_mpi 
       	= dynamic_cast<GroupProcessMpi*> (group_process);

      if (group_process_mpi != NULL) {
	group_process_mpi->set_send_blocking(blocking_send);
	group_process_mpi->set_recv_blocking(blocking_recv);
      }

#endif

      int rank_source = (rank+1)%size;
      int rank_dest   = (rank-1+size)%size;
      int array_size  = n*sizeof(double);

      void * handle_send = 
	group_process->send_begin (rank_source, array_source, array_size);
      void * handle_recv = 
	group_process->recv_begin (rank_dest,   array_dest,   array_size);

      group_process->wait(handle_recv);
      group_process->wait(handle_send);

      group_process->send_end(handle_send);
      group_process->recv_end(handle_recv);

      unit_assert(test_array(array_source,n+1,rank,rank));
      unit_assert(test_array(array_dest,  n+1,rank,rank_dest));

    }
  }

  for (int blocking_send = 0; blocking_send < num_blocking; blocking_send++) {
    for (int blocking_recv = 0; blocking_recv < num_blocking; blocking_recv++) {
      // MPI_Isend MPI_Irecv
      // MPI_Isend MPI_Recv
      // MPI_Send  MPI_Irecv
      // MPI_Send  MPI_Recv
      unit_func("GroupProcess","test");

      init_array(array_source,n+1,rank);
      init_array(array_dest,  n+1,rank);
      unit_assert(test_array(array_source,n+1,rank,rank));

#ifdef CONFIG_USE_MPI

      // Blocking send or receive is specific to MPI

      GroupProcessMpi * group_process_mpi
      	= dynamic_cast<GroupProcessMpi*> (group_process);

      group_process_mpi->set_send_blocking(blocking_send);
      group_process_mpi->set_recv_blocking(blocking_recv);

#endif

      int rank_source = (rank+1)%size;
      int rank_dest   = (rank-1+size)%size;
      int array_size  = n*sizeof(double);

      void * handle_send = 
	group_process->send_begin (rank_source, array_source, array_size);
      void * handle_recv = 
	group_process->recv_begin (rank_dest,   array_dest,   array_size);

      int counter = 0;
      while ( ! group_process->test(handle_recv) ||
	      ! group_process->test(handle_send) ) {
	// spinwait
	++ counter;
      }

      // assert recv completed and send completed

      group_process->send_end(handle_send);
      group_process->recv_end(handle_recv);

      unit_assert(test_array(array_source,n+1,rank,rank));
      unit_assert(test_array(array_dest,  n+1,rank,rank_dest));

    }
  }

  unit_func("GroupProcess","sync");
  switch (rank) {
  case 0:
    group_process->sync(1); // 0 - 1
    group_process->sync(2); // 0 - 2
    group_process->sync(3); // 0 - 3
    break;
  case 1:
    group_process->sync(0); // 0 - 1
    group_process->sync(3); // 1 - 3
    group_process->sync(2); // 1 - 2
    break;
  case 2:
    group_process->sync(3); // 2 - 3
    group_process->sync(0); // 0 - 2
    group_process->sync(1); // 1 - 2
    break;
  case 3:
    group_process->sync(2); // 2 - 3
    group_process->sync(1); // 1 - 3
    group_process->sync(0); // 0 - 3
    break;
  }
  unit_assert(true);

  // --------------------------------------------------
  unit_func("GroupProcess","barrier");
  group_process->barrier();
  delete group_process;
  
  unit_func("GroupProcess","bulk_send_add");
  unit_assert(unit_incomplete);
  unit_func("GroupProcess","bulk_send");
  unit_assert(unit_incomplete);
  unit_func("GroupProcess","bulk_send_wait");
  unit_assert(unit_incomplete);
  unit_func("GroupProcess","bulk_recv_add");
  unit_assert(unit_incomplete);
  unit_func("GroupProcess","bulk_recv");
  unit_assert(unit_incomplete);
  unit_func("GroupProcess","bulk_recv_wait");
  unit_assert(unit_incomplete);

  fflush(stdout);

  unit_finalize();

  PARALLEL_EXIT;

}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_GroupProcess.def.h)
