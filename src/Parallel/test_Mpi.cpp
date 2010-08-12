// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     test_Mpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 11:29:22 PDT 2010
/// @brief    Program implementing unit tests for the Mpi class


#include <stdio.h>
#include <string>

#include "cello.hpp"
#include "parallel.hpp"
#include "test.hpp"

int main(int argc, char ** argv)
{

  Mpi::init(&argc,&argv);

  unit_init (Mpi::rank(), Mpi::size());
  unit_class ("Mpi");
  unit_func("initialize");
  unit_assert(true);

  unit_func("size");
  int np = Mpi::size();
  unit_assert(np > 1);

  unit_func("rank");
  int ip = Mpi::rank();
  unit_assert(ip < np);
  unit_assert(ip >= 0);

  unit_func("is_root");
  unit_assert(Mpi::is_root() == (Mpi::rank()==0));

  unit_func("finalize");
  Mpi::finalize();
  unit_assert(true);

  unit_finalize();
}
