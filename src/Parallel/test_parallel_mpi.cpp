// $Id: test_parallel_mpi.cpp 1302 2010-03-17 00:16:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_parallel_mpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 11:29:22 PDT 2010
/// @brief    Program implementing unit tests for the ParallelMpi class
 
#include <stdio.h>
#include <string>

#include "parallel.hpp"
#include "test.hpp"


int main(int argc, char ** argv)
{
  unit_class ("ParallelMpi");

  unit_func("ParallelMpi");

  ParallelMpi mpi;
  unit_assert(true);

  unit_func("initialize");
  mpi.initialize(&argc,&argv);
  unit_assert(true);

  unit_func("size");
  printf ("size = %d\n",mpi.size());
  unit_assert(true);

  unit_func("rank");
  printf ("rank = %d\n",mpi.rank());
  unit_assert(true);

  unit_func("set_send_blocking");
  mpi.set_send_blocking(true);
  unit_assert (mpi.get_send_blocking() == true);
  unit_func("get_send_blocking");
  mpi.set_send_blocking(false);
  unit_assert (mpi.get_send_blocking() == false);

  unit_func("set_recv_blocking");
  mpi.set_recv_blocking(true);
  unit_assert (mpi.get_recv_blocking() == true);
  unit_func("get_recv_blocking");
  mpi.set_recv_blocking(false);
  unit_assert (mpi.get_recv_blocking() == false);

  unit_func("set_send_type");
  mpi.set_send_type(send_type_standard);
  unit_assert (mpi.get_send_type() == send_type_standard);
  unit_func("get_send_type");
  mpi.set_send_type(send_type_buffered);
  unit_assert (mpi.get_send_type() == send_type_buffered);
  mpi.set_send_type(send_type_ready);
  unit_assert (mpi.get_send_type() == send_type_ready);

  unit_func("finalize");
  mpi.finalize();
  unit_assert(true);
}
