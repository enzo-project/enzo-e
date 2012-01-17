// See LICENSE_CELLO file for license and copyright information

/// @file     test_Tree.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-07
/// @brief    Test program for the Tree class

#include "main.hpp"
#include "test.hpp"
#include "test_mesh.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("Tree");

  Tree * tree = test_tree_22();

  //--------------------------------------------------

  unit_func("Tree()");
  unit_assert (tree != NULL);

  //--------------------------------------------------
  unit_func("dimension()");
  unit_assert (tree->dimension() == 2);
  
  //--------------------------------------------------
  unit_func("refinement()");
  unit_assert (tree->refinement() == 2);

  //--------------------------------------------------
  unit_func("num_children()");
  unit_assert (tree->num_children() == 4);

  //--------------------------------------------------
  unit_func("num_nodes ()");
  unit_assert (tree->num_nodes() == 33);

  //--------------------------------------------------
  unit_func("root_node()");
  unit_assert (tree->root_node() != 0);

  //--------------------------------------------------
  unit_func("refine_node () (tested in test_ItNode)");
  unit_assert (true);

  //--------------------------------------------------
  unit_func("delete_node ()");

  unit_assert(tree->root_node()->child(1)->is_leaf() == false);

  const int max_depth = 4;
  NodeTrace node_trace (tree->root_node(), max_depth);
  
  node_trace.push(1);

  tree->delete_node(&node_trace);

  unit_assert(tree->num_nodes()==29);
  unit_assert(tree->root_node()->child(1)->is_leaf() == true);

  unit_assert(tree->root_node()->child(0)->is_leaf() == false);

  node_trace.pop();
  node_trace.push(0);
  tree->delete_node(&node_trace);

  unit_assert(tree->num_nodes()==17);
  unit_assert(tree->root_node()->child(0)->is_leaf() == true);
  
  //--------------------------------------------------

  unit_func("balance_node ()");
  unit_assert (false);

  //--------------------------------------------------
  unit_func("node_neighbor ()");
  unit_assert (false);

  //--------------------------------------------------
  unit_func("node_parent ()");
  unit_assert (false);

  //--------------------------------------------------
  unit_func("node_child ()");
  unit_assert (false);

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

