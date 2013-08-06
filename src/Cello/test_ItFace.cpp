// See LICENSE_CELLO file for license and copyright information

/// @file     test_ItFace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the ItFace class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ItFace");

  int count;
  int if3[3];

  count = 0;
  unit_func("ItFace(1,0)");
  ItFace it_face10 (1,0); 
  while (it_face10.next(if3)) {
    ++count;
  }
  unit_assert(count == 3 - 1);

  count = 0;
  unit_func("ItFace(2,0)");
  ItFace it_face20 (2,0); 
  while (it_face20.next(if3)) {
    ++count;
  }
  unit_assert(count == 9 - 1);

  count = 0;
  unit_func("ItFace(2,1)");
  ItFace it_face21 (2,1); 
  while (it_face21.next(if3)) {
    ++count;
  }
  unit_assert(count == 9 - 4 - 1);

  count = 0;
  unit_func("ItFace(3,0)");
  ItFace it_face30 (3,0); 
  while (it_face30.next(if3)) {
    ++count;
  }
  unit_assert(count == 27 - 1);

  count = 0;
  unit_func("ItFace(3,1)");
  ItFace it_face31 (3,1); 
  while (it_face31.next(if3)) {
    ++count;
  }
  unit_assert(count == 27 - 8 - 1);

  count = 0;
  unit_func("ItFace(3,2)");
  ItFace it_face32 (3,2); 
  while (it_face32.next(if3)) {
    ++count;
  }
  unit_assert(count == 27 - 8 - 12 - 1);

  //--------------------------------------------------

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

