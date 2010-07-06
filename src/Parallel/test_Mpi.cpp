// $Id: test_parallel_mpi.cpp 1302 2010-03-17 00:16:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 11:29:22 PDT 2010
/// @brief    Program implementing unit tests for the ParallelMpi class

#include <mpi.h> 
#include <stdio.h>
#include <string>

#include "cello.hpp"
#include "parallel.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{
  ParallelMpi * parallel = new ParallelMpi;

  parallel->initialize(&argc,&argv);
  unit_init (parallel->process_rank(), parallel->process_count());
  unit_func("initialize");
  unit_class ("ParallelMpi");
  unit_assert(true);

  unit_func("process_count");
  int np = parallel->process_count();
  unit_assert(np > 1);

  unit_func("process_rank");
  int ip = parallel->process_rank();
  unit_assert(ip < np);
  unit_assert(ip >= 0);

  unit_func("set_send_blocking");
  parallel->set_send_blocking(true);
  unit_assert (parallel->get_send_blocking() == true);
  unit_func("get_send_blocking");
  parallel->set_send_blocking(false);
  unit_assert (parallel->get_send_blocking() == false);

  unit_func("set_recv_blocking");
  parallel->set_recv_blocking(true);
  unit_assert (parallel->get_recv_blocking() == true);
  unit_func("get_recv_blocking");
  parallel->set_recv_blocking(false);
  unit_assert (parallel->get_recv_blocking() == false);

  unit_func("set_send_type");
  parallel->set_send_type(send_standard);
  unit_assert (parallel->get_send_type() == send_standard);
  unit_func("get_send_type");
  parallel->set_send_type(send_buffered);
  unit_assert (parallel->get_send_type() == send_buffered);
  parallel->set_send_type(send_ready);
  unit_assert (parallel->get_send_type() == send_ready);

  unit_func("finalize");
  parallel->finalize();
  unit_assert(true);

  unit_finalize();
}
