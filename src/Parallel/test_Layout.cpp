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

#define INDEX3(I,N) I[0] + N[0]*(I[1] + N[1]*I[2])

int main(int argc, char ** argv)
{
  unit_class("Layout");
  unit_func("Layout");
  Layout layout_serial();
  unit_assert (true);
}
