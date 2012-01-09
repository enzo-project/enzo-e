// See LICENSE_CELLO file for license and copyright information

/// @file     test_Tree.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the Tree class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Tree");

  int d = 2;
  int r = 2;

  Tree * tree = new Tree (d,r);

  unit_func("Tree");
  unit_assert (tree != NULL);

  unit_func("dimension");
  unit_assert (false);
  
  unit_func("refinement");
  unit_assert (false);

  unit_func("num_nodes ");
  unit_assert (false);

  unit_func("num_levels ");
  unit_assert (false);

  unit_func("root_node");
  unit_assert (false);

  unit_func("find_node ");
  unit_assert (false);

  unit_func("refine_node ");
  unit_assert (false);

  unit_func("delete_node ");
  unit_assert (false);

  unit_func("balance_node ");
  unit_assert (false);

  unit_func("node_neighbor ");
  unit_assert (false);

  unit_func("node_parent ");
  unit_assert (false);

  unit_func("node_child ");
  unit_assert (false);

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

