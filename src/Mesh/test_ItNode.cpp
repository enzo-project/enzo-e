// See LICENSE_CELLO file for license and copyright information

/// @file     test_ItNode.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-12
/// @brief    Test program for the ItNode class

#include "main.hpp"
#include "test.hpp"
#include "test_mesh.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ItNode");

  // Setup Tree

  Tree * tree = test_tree_22();
  int max_depth = 3;

  // Test ItNode depth-first Morton ordering

  ItNode it_node(tree, max_depth);

  Node * root = tree->root_node();

  for (int count = 0; count < 2; count ++) {
    unit_assert (++it_node == root);
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(0));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(1));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(1)->child(0));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(1)->child(1));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(1)->child(2));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(1)->child(3));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(2));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(2)->child(0));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(2)->child(1));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(2)->child(2));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(2)->child(3));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(0)->child(3));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(1));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(1)->child(0));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(1)->child(1));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(1)->child(2));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(1)->child(3));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(2));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0)->child(0));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0)->child(1));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0)->child(2));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0)->child(3));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(1));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(2));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(3));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0)->child(0));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0)->child(1));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0)->child(2));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == root->child(3)->child(0)->child(3));
    unit_assert (it_node.done()==false);
    unit_assert (++it_node == 0);
    unit_assert (it_node.done()==true);
  }

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

