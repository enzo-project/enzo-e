// See LICENSE_CELLO file for license and copyright information

/// @file     test_NodeInfo.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the NodeInfo class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("NodeInfo");

  int d = 3;
  int r = 2;
  //  NodeInfo * NodeInfo = new NodeInfo (d,r);
  NodeInfo * node_info = new NodeInfo (d,r);

  unit_assert (node_info != NULL);

  //--------------------------------------------------

  unit_func ("trace()");

  unit_assert (false);

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

