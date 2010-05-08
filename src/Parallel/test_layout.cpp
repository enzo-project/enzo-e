// $Id: test_block.cpp 1369 2010-04-08 01:38:06Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-19
/// @brief    Unit tests for the Layout class
///
/// Run with mpirun -np 4


#include <mpi.h>

#include "cello.h"

#include "error.hpp"
#include "test.hpp"
#include "parallel.hpp"
#include "array.hpp"

#define INDEX3(I,N) I[0] + N[0]*(I[1] + N[1]*I[2])

int main(int argc, char ** argv)
{

  Parallel * parallel = Parallel::instance();

  parallel->initialize(&argc, &argv);

  int process_count = parallel->process_count();
  int process_rank  = parallel->process_rank();

  if (process_count != 4) {
    printf ("mpirun -np 4 required\n");
    parallel->finalize();
    parallel->abort();
  }

  // unit_note("mpirun -np 8 required");
  // unit_assert(process_count == 8);

  // Initialize unit testing for parallel runs

  unit_init (process_rank, process_count);

  // Set unit testing class

  unit_class ("Layout");

  unit_size<Layout> ();

  // Initialize test-independent problem parameters

  int dim = 3;             // Problem dimension dim = 3
  int n3[]   = {15,7,19};    // Array size n3[]
  int np3[3] = { 5,2, 3};    // Process block size
  int process_block_count = np3[0]*np3[1]*np3[2];

  //----------------------------------------------------------------------
  // TEST 1:  One (P=1,T=1) [Serial]
  //----------------------------------------------------------------------

  // create serial layouts on each process
  unit_func("Layout");
  Layout * layout_serial = new Layout (dim);
  unit_assert (true);

  // set the layout's array size and confirm
  unit_func("set_array");
  layout_serial -> set_array (dim, n3);

  int n3_serial[3];
  unit_func("array_size");
  layout_serial->array_size(3,n3_serial);
  unit_assert (n3[0]==n3_serial[0] &&
	       n3[1]==n3_serial[1] &&
	       n3[2]==n3_serial[2]);

  // set the layout's process range (process_rank) and confirm

  unit_func("set_processes");
  layout_serial->set_processes(process_rank,1);

  unit_func("process_first");
  unit_assert (layout_serial -> process_first() == process_rank);
  unit_func("process_count");
  unit_assert (layout_serial -> process_count() == 1);

  // set the layout's thread range and confirm

  unit_func("set_threads");
  layout_serial->set_threads(0,1);   

  unit_func("thread_first");
  unit_assert (layout_serial -> thread_first() == 0);
  unit_func("thread_count");
  unit_assert (layout_serial -> thread_count() == 1);

  // set the layout's process block size

  layout_serial -> set_process_blocks (3,np3);

  // confirm number of process blocks assigned to this process

  unit_func("process_block_count");
  int expected_serial = process_block_count;

  unit_assert (layout_serial -> process_block_count(process_rank) == 
	       expected_serial);
  
  // confirm indices of process blocks assigned to this process

  int test_index;
  int computed_index;
  int process_index[3];

  bool success = true;
  for (test_index = 0; 
       test_index < expected_serial; 
       test_index ++) 
    {

      // Get process_index[] array
      layout_serial ->process_block_indices(3, process_rank,test_index, process_index);

      // back-compute index from array
      computed_index = INDEX3(process_index,np3);

      // ensure that back-computed index is the same as the test index
      if (! computed_index == test_index) {
	unit_assert (success=false);
      }

      // compare with process_block_number(), the inverse operation of
      // process_block_indices()
      int ip=-1,ipb=-1;
      layout_serial -> process_block_number (dim, &ip,&ipb, process_index);
      if (! ipb == test_index) {
	unit_assert (success=false);
      }
    }

  unit_assert (success);

  unit_func("~Layout");
  delete layout_serial;
  unit_assert (true);
  
  //----------------------------------------------------------------------
  // TEST 2:  One (P=p,T=1) [MPI parallel]
  //----------------------------------------------------------------------

  // create mpi parallel layouts on each process
  unit_func("Layout");
  Layout * layout_mpi = new Layout (dim);
  unit_assert (true);

  // set the layout's array size and confirm
  unit_func("set_array");
  layout_mpi -> set_array (dim, n3);

  int n3_mpi[3];
  unit_func("array_size");
  layout_mpi->array_size(3,n3_mpi);
  unit_assert (n3[0]==n3_mpi[0] &&
	       n3[1]==n3_mpi[1] &&
	       n3[2]==n3_mpi[2]);

  // set the layout's process range (process_rank) and confirm

  unit_func("set_processes");
  layout_mpi->set_processes(0,process_count);

  unit_func("process_first");
  unit_assert (layout_mpi -> process_first() == 0);

  unit_func("process_count");
  unit_assert (layout_mpi -> process_count() == process_count);

  // set the layout's thread range and confirm

  unit_func("set_threads");
  layout_mpi->set_threads(0,1);   

  unit_func("thread_first");
  unit_assert (layout_mpi -> thread_first() == 0);
  unit_func("thread_count");
  unit_assert (layout_mpi -> thread_count() == 1);

  // set the layout's process block size

  layout_mpi -> set_process_blocks (3,np3);

  // confirm number of process blocks assigned to this process

  unit_func("process_block_count");
  int expected_mpi = 
    (process_rank + 1) * process_block_count / process_count -
    (process_rank + 0) * process_block_count / process_count;
    
  unit_assert (layout_mpi -> process_block_count(process_rank) == 
	       expected_mpi);
  
  // confirm indices of process blocks assigned to this process

  
  success = true;

  int offset_mpi = process_rank * process_block_count / process_count;

  for (test_index = 0; 
       test_index < expected_mpi; 
       test_index ++) 
    {

      // Get process_index[] array
      layout_mpi ->process_block_indices(3, process_rank,test_index, process_index);

      // back-compute index from array
      computed_index = INDEX3(process_index,np3);

      // ensure that back-computed index is the same as the test index
      if (! computed_index == test_index + offset_mpi) {
	unit_assert (success=false);
      }

      // compare with process_block_number(), the inverse operation of
      // process_block_indices()
      int ip=-1,ipb=-1;
      layout_mpi -> process_block_number (dim, &ip,&ipb, process_index);
      if (! ipb == test_index) {
	unit_assert (success=false);
      }
    }

  unit_assert (success);

  unit_func("~Layout");
  delete layout_mpi;
  unit_assert (true);

  //----------------------------------------------------------------------
  // TEST 3:  Two (P=p,T=t) [Multiple layouts]
  //----------------------------------------------------------------------

  parallel->finalize();

}
