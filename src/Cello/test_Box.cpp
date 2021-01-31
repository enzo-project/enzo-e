// See LICENSE_CELLO file for license and copyright information

/// @file     test_Box.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2021-01-28
/// @brief    Test program for the Box class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Box");

  int n3[] = {4,8,6};
  int g3[] = {2,1,0};
  Box * box = new Box (3,n3,g3);

  unit_assert (box != NULL);

  //--------------------------------------------------


  //--------------------------------------------------

  delete box;

  unit_finalize();

  exit_();
}

PARALLEL_MAIN_END

