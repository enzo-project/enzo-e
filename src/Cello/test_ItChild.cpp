// See LICENSE_CELLO file for license and copyright information

/// @file     test_ItChild.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the ItChild class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ItChild");

  int count;
  int ic3[3];

  count = 0;
  unit_func("ItChild(1)");
  ItChild it_child10 (1); 
  while (it_child10.next(ic3)) {
    ++count;
  }
  unit_assert(count == 2);

  count = 0;
  unit_func("ItChild(2)");
  ItChild it_child20 (2); 
  while (it_child20.next(ic3)) {
    ++count;
  }
  unit_assert(count == 4);

  count = 0;
  unit_func("ItChild(3)");
  ItChild it_child21 (3); 
  while (it_child21.next(ic3)) {
    ++count;
  }
  unit_assert(count == 8);

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

