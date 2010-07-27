// $Id: test_Mpi.cpp 1302 2010-03-17 00:16:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 11:29:22 PDT 2010
/// @brief    Program implementing unit tests for the Mpi class

#include <mpi.h> 
#include <stdio.h>
#include <string>

#include "cello.hpp"
#include "parallel.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{
  GroupProcessMpi process_group;
  process_group.initialize(&argc,&argv);
  unit_init (process_group.rank(), process_group.size());
  unit_class ("Mpi");
  unit_func("initialize");
  unit_assert(true);

  unit_func("size");
  int np = process_group.size();
  unit_assert(np > 1);

  unit_func("rank");
  int ip = process_group.rank();
  unit_assert(ip < np);
  unit_assert(ip >= 0);

  unit_func("is_root");
  unit_assert(process_group.is_root() == (process_group.rank()==0));

  unit_func("finalize");
  process_group.finalize();
  unit_assert(true);

  unit_finalize();
}
