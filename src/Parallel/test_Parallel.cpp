// $Id: test_Parallel.cpp 1372 2010-04-08 05:36:42Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Parallel.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @bug      Crashes in Parallel::initialize() in MPI_Init with LAM MPI
/// @date     Tue Apr 20 14:19:04 PDT 2010
/// @brief    Program implementing unit tests for the Parallel class

#include "mpi.h" 
#include "cello.hpp"
#include "parallel.hpp"
#include "test.hpp"

void init_array(double * array, int length, int rank)
{
  for (int i=0; i<length-1; i++) {
    array[i] = 1.0*rank;
  }
  array[length-1] = -1.0*rank;
}

bool test_array(double * array, int length, int rank, int value)
{
  bool valid = true;
  for (int i=0; i<length-1; i++) {
    valid = valid && (array[i] == 1.0*value);
  }
  valid = valid && (array[length-1] == -1.0*rank);
  if (!valid) {
    printf ("%d expected [%g %g %g %g %g]; actual [%g %g %g %g %g]\n",
	    rank,
	    1.0*value,1.0*value,1.0*value,1.0*value,-1.0*rank,
	    array[0], array[1], array[2], array[3], array[4]);
  }
  return valid;
}

int main(int argc, char ** argv)
{
  MPI_Init(&argc,&argv);

  //  int buffer[500];
  //  for (int i=0; i<500; i++) buffer[i] = 0; 

  GroupProcessMpi process_group;


  int rank = process_group.rank();
  int size = process_group.size();

  unit_init(rank,size);

  unit_class("GroupProcessMpi");

  unit_func("size");
  unit_assert(size == 4);

  unit_func("rank");
  unit_assert(0 <= rank && rank < 4);

  unit_func("is_root");
  unit_assert(process_group.is_root() == (rank == 0));

  // Test that init_array() and test_array() work independently of Parallel
  const int n = 4;
  double array_source[n+1], array_dest[n+1];
  init_array(array_source,n+1,rank);
  init_array(array_dest,  n+1,rank);
  unit_assert(test_array(array_source,n+1,rank,rank));

  unit_func("send");

  int rank_source = (rank+1)%size;
  int rank_dest   = (rank-1+size)%size;
  int array_size  = n*sizeof(double);

  printf ("%d %d %d %d\n",rank,size,rank_source,rank_dest);

  int handle_send = process_group.send(rank_source, array_source, array_size);
  int handle_recv = process_group.recv(rank_dest,   array_dest,   array_size);

  //  for (int i=0; i<500; i++) if (buffer[i] != 0) printf ("X: %d %d\n",i,buffer[i]); 

  MPI_Finalize();
  exit(0);
  //  int i,j;

//   process_group.recv_wait(handle_recv);
//  process_group.send_wait(handle_send);

//   unit_assert(test_array(array_source,n+1,rank,rank));
//   unit_assert(test_array(array_dest,  n+1,rank,rank_dest));

//   unit_func("barrier");
//   process_group.barrier();
  
//   unit_func("wait");
//   switch (rank) {
//   case 0:
//     process_group.wait(1); // 0 - 1
//     process_group.wait(2); // 0 - 2
//     process_group.wait(3); // 0 - 3
//     break;
//   case 1:
//     process_group.wait(0); // 0 - 1
//     process_group.wait(3); // 1 - 3
//     process_group.wait(2); // 1 - 2
//     break;
//   case 2:
//     process_group.wait(3); // 2 - 3
//     process_group.wait(0); // 0 - 2
//     process_group.wait(1); // 1 - 2
//     break;
//   case 3:
//     process_group.wait(2); // 2 - 3
//     process_group.wait(1); // 1 - 3
//     process_group.wait(0); // 0 - 3
//     break;
//   }
//   unit_assert(true);

//   unit_func("bulk_send_add");
//   unit_assert(false);
//   unit_func("bulk_send");
//   unit_assert(false);
//   unit_func("bulk_send_wait");
//   unit_assert(false);
//   unit_func("bulk_recv_add");
//   unit_assert(false);
//   unit_func("bulk_recv");
//   unit_assert(false);
//   unit_func("bulk_recv_wait");
//   unit_assert(false);
  
//   unit_finalize();
//   Mpi::finalize();

}
