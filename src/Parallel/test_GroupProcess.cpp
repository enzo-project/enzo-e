// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_GroupProcess.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      Crashes in Parallel::initialize() in MPI_Init with LAM MPI
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
    char message[ERROR_MESSAGE_LENGTH];
    sprintf (message,
	     "%d expected [%g %g %g %g %g]; actual [%g %g %g %g %g]\n",
	     rank,
	     1.0*value,1.0*value,1.0*value,1.0*value,-1.0*rank,
	     array[0], array[1], array[2], array[3], array[4]);
    WARNING_MESSAGE("test_array",message);
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


  GroupProcess * parallel = GroupProcess::create();

  int rank = parallel->rank();
  int size = parallel->size();

  unit_init(parallel->rank(), parallel->size());

  unit_class("GroupProcess");

  unit_func("size");
  unit_assert(size == np);

  unit_func("rank");
  unit_assert(0 <= rank && rank < np);

  unit_func("is_root");
  unit_assert(parallel->is_root() == (rank == 0));

  // Test that init_array() and test_array() work independently of Parallel
  const int n = np;
  double array_source[n+1], array_dest[n+1];

  unit_func("wait");

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

      GroupProcessMpi * parallel_mpi 
       	= dynamic_cast<GroupProcessMpi*> (parallel);

      parallel_mpi->set_send_blocking(blocking_send);
      parallel_mpi->set_recv_blocking(blocking_recv);

#endif

      int rank_source = (rank+1)%size;
      int rank_dest   = (rank-1+size)%size;
      int array_size  = n*sizeof(double);

      void * handle_send = 
	parallel->send_begin (rank_source, array_source, array_size);
      void * handle_recv = 
	parallel->recv_begin (rank_dest,   array_dest,   array_size);

      parallel->wait(handle_recv);
      parallel->wait(handle_send);

      parallel->send_end(handle_send);
      parallel->recv_end(handle_recv);

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
      unit_func("test");

      init_array(array_source,n+1,rank);
      init_array(array_dest,  n+1,rank);
      unit_assert(test_array(array_source,n+1,rank,rank));

#ifdef CONFIG_USE_MPI

      // Blocking send or receive is specific to MPI

      GroupProcessMpi * parallel_mpi
      	= dynamic_cast<GroupProcessMpi*> (parallel);

      parallel_mpi->set_send_blocking(blocking_send);
      parallel_mpi->set_recv_blocking(blocking_recv);

#endif

      int rank_source = (rank+1)%size;
      int rank_dest   = (rank-1+size)%size;
      int array_size  = n*sizeof(double);

      void * handle_send = 
	parallel->send_begin (rank_source, array_source, array_size);
      void * handle_recv = 
	parallel->recv_begin (rank_dest,   array_dest,   array_size);

      int counter = 0;
      while ( ! parallel->test(handle_recv) ||
	      ! parallel->test(handle_send) ) {
	// spinwait
	++ counter;
      }

      // assert recv completed and send completed

      parallel->send_end(handle_send);
      parallel->recv_end(handle_recv);

      unit_func("test");

      unit_assert(test_array(array_source,n+1,rank,rank));
      unit_assert(test_array(array_dest,  n+1,rank,rank_dest));

    }
  }

  unit_func("sync");
  switch (rank) {
  case 0:
    parallel->sync(1); // 0 - 1
    parallel->sync(2); // 0 - 2
    parallel->sync(3); // 0 - 3
    break;
  case 1:
    parallel->sync(0); // 0 - 1
    parallel->sync(3); // 1 - 3
    parallel->sync(2); // 1 - 2
    break;
  case 2:
    parallel->sync(3); // 2 - 3
    parallel->sync(0); // 0 - 2
    parallel->sync(1); // 1 - 2
    break;
  case 3:
    parallel->sync(2); // 2 - 3
    parallel->sync(1); // 1 - 3
    parallel->sync(0); // 0 - 3
    break;
  }
  unit_assert(true);

  // --------------------------------------------------
  unit_func("barrier");
  parallel->barrier();
  delete parallel;
  
  unit_func("bulk_send_add");
  unit_assert(false);
  unit_func("bulk_send");
  unit_assert(false);
  unit_func("bulk_send_wait");
  unit_assert(false);
  unit_func("bulk_recv_add");
  unit_assert(false);
  unit_func("bulk_recv");
  unit_assert(false);
  unit_func("bulk_recv_wait");
  unit_assert(false);

  fflush(stdout);

  unit_finalize();

  PARALLEL_EXIT;

}

PARALLEL_MAIN_END

#include PARALLEL_CHARM_INCLUDE(test_GroupProcess.def.h)
