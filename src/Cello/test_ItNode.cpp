// See LICENSE_CELLO file for license and copyright information

/// @file     test_ItNode.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-12
/// @brief    Test program for the ItNode class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"
#include "mesh_functions.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("ItNode");

  //--------------------------------------------------
  // Tree (2,2)
  //--------------------------------------------------

  {
    Tree * tree = test_tree_22();

    tree_to_png (tree, "test_ItNode.png",512,512);

    // Test ItNode depth-first Morton ordering

    ItNode it_node(tree);

    Node * root = tree->root_node();

    for (int count = 0; count < 2; count ++) {
      unit_assert (it_node.next_leaf() == root->child(0)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(1)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(1)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(1)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(1)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(2)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(2)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(2)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(2)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(0)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(1)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(1)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(1)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(1)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(0)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(0)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(0)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(0)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(3)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(3)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(3)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == root->child(3)->child(3)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_leaf() == 0);
      unit_assert (it_node.done()==true);
    }

    for (int count = 0; count < 2; count ++) {
      unit_assert (it_node.next_node() == root);
      printf ("%p\n",it_node.node_trace()->node());
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0));
      printf ("%p\n",it_node.node_trace()->node());
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(0));
      printf ("%p\n",it_node.node_trace()->node());
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(1)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(1)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(1)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(1)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(2)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(2)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(2)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(2)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(0)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(1)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(1)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(1)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(1)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(0)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(0)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(0)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(0)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(3)->child(0));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(3)->child(1));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(3)->child(2));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == root->child(3)->child(3)->child(3));
      unit_assert (it_node.done()==false);
      unit_assert (it_node.next_node() == 0);
      unit_assert (it_node.done()==true);
    }
  }

  //--------------------------------------------------
  // Tree (3,2)
  //--------------------------------------------------

  {

    Tree * tree = test_tree_32();

    // Test ItNode depth-first Morton ordering

    ItNode it_node(tree);

    Node * root = tree->root_node();

    for (int count = 0; count < 2; count ++) {
      unit_assert (it_node.next_leaf() == root->child(0)->child(0));
      unit_assert (it_node.next_leaf() == root->child(0)->child(1));
      unit_assert (it_node.next_leaf() == root->child(0)->child(2));
      unit_assert (it_node.next_leaf() == root->child(0)->child(3));
      unit_assert (it_node.next_leaf() == root->child(0)->child(4));
      unit_assert (it_node.next_leaf() == root->child(0)->child(5));
      unit_assert (it_node.next_leaf() == root->child(0)->child(6));
      unit_assert (it_node.next_leaf() == root->child(0)->child(7)->child(0));
      unit_assert (it_node.next_leaf() == root->child(0)->child(7)->child(1));
      unit_assert (it_node.next_leaf() == root->child(0)->child(7)->child(2));
      unit_assert (it_node.next_leaf() == root->child(0)->child(7)->child(3));
      unit_assert (it_node.next_leaf() == root->child(0)->child(7)->child(4));
      unit_assert (it_node.next_leaf() == root->child(0)->child(7)->child(5));
      unit_assert (it_node.next_leaf() == root->child(0)->child(7)->child(6));
      unit_assert (it_node.next_leaf() == root->child(0)->child(7)->child(7));
      unit_assert (it_node.next_leaf() == root->child(1)->child(0));
      unit_assert (it_node.next_leaf() == root->child(1)->child(1));
      unit_assert (it_node.next_leaf() == root->child(1)->child(2));
      unit_assert (it_node.next_leaf() == root->child(1)->child(3));
      unit_assert (it_node.next_leaf() == root->child(1)->child(4));
      unit_assert (it_node.next_leaf() == root->child(1)->child(5));
      unit_assert (it_node.next_leaf() == root->child(1)->child(6)->child(0));
      unit_assert (it_node.next_leaf() == root->child(1)->child(6)->child(1));
      unit_assert (it_node.next_leaf() == root->child(1)->child(6)->child(2));
      unit_assert (it_node.next_leaf() == root->child(1)->child(6)->child(3));
      unit_assert (it_node.next_leaf() == root->child(1)->child(6)->child(4));
      unit_assert (it_node.next_leaf() == root->child(1)->child(6)->child(5));
      unit_assert (it_node.next_leaf() == root->child(1)->child(6)->child(6));
      unit_assert (it_node.next_leaf() == root->child(1)->child(6)->child(7));
      unit_assert (it_node.next_leaf() == root->child(1)->child(7));
      unit_assert (it_node.next_leaf() == root->child(2)->child(0));
      unit_assert (it_node.next_leaf() == root->child(2)->child(1));
      unit_assert (it_node.next_leaf() == root->child(2)->child(2));
      unit_assert (it_node.next_leaf() == root->child(2)->child(3));
      unit_assert (it_node.next_leaf() == root->child(2)->child(4));
      unit_assert (it_node.next_leaf() == root->child(2)->child(5)->child(0));
      unit_assert (it_node.next_leaf() == root->child(2)->child(5)->child(1));
      unit_assert (it_node.next_leaf() == root->child(2)->child(5)->child(2));
      unit_assert (it_node.next_leaf() == root->child(2)->child(5)->child(3));
      unit_assert (it_node.next_leaf() == root->child(2)->child(5)->child(4));
      unit_assert (it_node.next_leaf() == root->child(2)->child(5)->child(5));
      unit_assert (it_node.next_leaf() == root->child(2)->child(5)->child(6));
      unit_assert (it_node.next_leaf() == root->child(2)->child(5)->child(7));
      unit_assert (it_node.next_leaf() == root->child(2)->child(6));
      unit_assert (it_node.next_leaf() == root->child(2)->child(7));
      unit_assert (it_node.next_leaf() == root->child(3)->child(0));
      unit_assert (it_node.next_leaf() == root->child(3)->child(1));
      unit_assert (it_node.next_leaf() == root->child(3)->child(2));
      unit_assert (it_node.next_leaf() == root->child(3)->child(3));
      unit_assert (it_node.next_leaf() == root->child(3)->child(4)->child(0));
      unit_assert (it_node.next_leaf() == root->child(3)->child(4)->child(1));
      unit_assert (it_node.next_leaf() == root->child(3)->child(4)->child(2));
      unit_assert (it_node.next_leaf() == root->child(3)->child(4)->child(3));
      unit_assert (it_node.next_leaf() == root->child(3)->child(4)->child(4));
      unit_assert (it_node.next_leaf() == root->child(3)->child(4)->child(5));
      unit_assert (it_node.next_leaf() == root->child(3)->child(4)->child(6));
      unit_assert (it_node.next_leaf() == root->child(3)->child(4)->child(7));
      unit_assert (it_node.next_leaf() == root->child(3)->child(5));
      unit_assert (it_node.next_leaf() == root->child(3)->child(6));
      unit_assert (it_node.next_leaf() == root->child(3)->child(7));
      unit_assert (it_node.next_leaf() == root->child(4)->child(0));
      unit_assert (it_node.next_leaf() == root->child(4)->child(1));
      unit_assert (it_node.next_leaf() == root->child(4)->child(2));
      unit_assert (it_node.next_leaf() == root->child(4)->child(3)->child(0));
      unit_assert (it_node.next_leaf() == root->child(4)->child(3)->child(1));
      unit_assert (it_node.next_leaf() == root->child(4)->child(3)->child(2));
      unit_assert (it_node.next_leaf() == root->child(4)->child(3)->child(3));
      unit_assert (it_node.next_leaf() == root->child(4)->child(3)->child(4));
      unit_assert (it_node.next_leaf() == root->child(4)->child(3)->child(5));
      unit_assert (it_node.next_leaf() == root->child(4)->child(3)->child(6));
      unit_assert (it_node.next_leaf() == root->child(4)->child(3)->child(7));
      unit_assert (it_node.next_leaf() == root->child(4)->child(4));
      unit_assert (it_node.next_leaf() == root->child(4)->child(5));
      unit_assert (it_node.next_leaf() == root->child(4)->child(6));
      unit_assert (it_node.next_leaf() == root->child(4)->child(7));
      unit_assert (it_node.next_leaf() == root->child(5)->child(0));
      unit_assert (it_node.next_leaf() == root->child(5)->child(1));
      unit_assert (it_node.next_leaf() == root->child(5)->child(2)->child(0));
      unit_assert (it_node.next_leaf() == root->child(5)->child(2)->child(1));
      unit_assert (it_node.next_leaf() == root->child(5)->child(2)->child(2));
      unit_assert (it_node.next_leaf() == root->child(5)->child(2)->child(3));
      unit_assert (it_node.next_leaf() == root->child(5)->child(2)->child(4));
      unit_assert (it_node.next_leaf() == root->child(5)->child(2)->child(5));
      unit_assert (it_node.next_leaf() == root->child(5)->child(2)->child(6));
      unit_assert (it_node.next_leaf() == root->child(5)->child(2)->child(7));
      unit_assert (it_node.next_leaf() == root->child(5)->child(3));
      unit_assert (it_node.next_leaf() == root->child(5)->child(4));
      unit_assert (it_node.next_leaf() == root->child(5)->child(5));
      unit_assert (it_node.next_leaf() == root->child(5)->child(6));
      unit_assert (it_node.next_leaf() == root->child(5)->child(7));
      unit_assert (it_node.next_leaf() == root->child(6)->child(0));
      unit_assert (it_node.next_leaf() == root->child(6)->child(1)->child(0));
      unit_assert (it_node.next_leaf() == root->child(6)->child(1)->child(1));
      unit_assert (it_node.next_leaf() == root->child(6)->child(1)->child(2));
      unit_assert (it_node.next_leaf() == root->child(6)->child(1)->child(3));
      unit_assert (it_node.next_leaf() == root->child(6)->child(1)->child(4));
      unit_assert (it_node.next_leaf() == root->child(6)->child(1)->child(5));
      unit_assert (it_node.next_leaf() == root->child(6)->child(1)->child(6));
      unit_assert (it_node.next_leaf() == root->child(6)->child(1)->child(7));
      unit_assert (it_node.next_leaf() == root->child(6)->child(2));
      unit_assert (it_node.next_leaf() == root->child(6)->child(3));
      unit_assert (it_node.next_leaf() == root->child(6)->child(4));
      unit_assert (it_node.next_leaf() == root->child(6)->child(5));
      unit_assert (it_node.next_leaf() == root->child(6)->child(6));
      unit_assert (it_node.next_leaf() == root->child(6)->child(7));
      unit_assert (it_node.next_leaf() == root->child(7)->child(0)->child(0));
      unit_assert (it_node.next_leaf() == root->child(7)->child(0)->child(1));
      unit_assert (it_node.next_leaf() == root->child(7)->child(0)->child(2));
      unit_assert (it_node.next_leaf() == root->child(7)->child(0)->child(3));
      unit_assert (it_node.next_leaf() == root->child(7)->child(0)->child(4));
      unit_assert (it_node.next_leaf() == root->child(7)->child(0)->child(5));
      unit_assert (it_node.next_leaf() == root->child(7)->child(0)->child(6));
      unit_assert (it_node.next_leaf() == root->child(7)->child(0)->child(7));
      unit_assert (it_node.next_leaf() == root->child(7)->child(1));
      unit_assert (it_node.next_leaf() == root->child(7)->child(2));
      unit_assert (it_node.next_leaf() == root->child(7)->child(3));
      unit_assert (it_node.next_leaf() == root->child(7)->child(4));
      unit_assert (it_node.next_leaf() == root->child(7)->child(5));
      unit_assert (it_node.next_leaf() == root->child(7)->child(6));
      unit_assert (it_node.next_leaf() == root->child(7)->child(7));
      unit_assert (it_node.next_leaf() == 0);
      unit_assert (it_node.done()==true);
    }

    for (int count = 0; count < 2; count ++) {
      unit_assert (it_node.next_node() == root);
      unit_assert (it_node.next_node() == root->child(0));
      unit_assert (it_node.next_node() == root->child(0)->child(0));
      unit_assert (it_node.next_node() == root->child(0)->child(1));
      unit_assert (it_node.next_node() == root->child(0)->child(2));
      unit_assert (it_node.next_node() == root->child(0)->child(3));
      unit_assert (it_node.next_node() == root->child(0)->child(4));
      unit_assert (it_node.next_node() == root->child(0)->child(5));
      unit_assert (it_node.next_node() == root->child(0)->child(6));
      unit_assert (it_node.next_node() == root->child(0)->child(7));
      unit_assert (it_node.next_node() == root->child(0)->child(7)->child(0));
      unit_assert (it_node.next_node() == root->child(0)->child(7)->child(1));
      unit_assert (it_node.next_node() == root->child(0)->child(7)->child(2));
      unit_assert (it_node.next_node() == root->child(0)->child(7)->child(3));
      unit_assert (it_node.next_node() == root->child(0)->child(7)->child(4));
      unit_assert (it_node.next_node() == root->child(0)->child(7)->child(5));
      unit_assert (it_node.next_node() == root->child(0)->child(7)->child(6));
      unit_assert (it_node.next_node() == root->child(0)->child(7)->child(7));
      unit_assert (it_node.next_node() == root->child(1));
      unit_assert (it_node.next_node() == root->child(1)->child(0));
      unit_assert (it_node.next_node() == root->child(1)->child(1));
      unit_assert (it_node.next_node() == root->child(1)->child(2));
      unit_assert (it_node.next_node() == root->child(1)->child(3));
      unit_assert (it_node.next_node() == root->child(1)->child(4));
      unit_assert (it_node.next_node() == root->child(1)->child(5));
      unit_assert (it_node.next_node() == root->child(1)->child(6));
      unit_assert (it_node.next_node() == root->child(1)->child(6)->child(0));
      unit_assert (it_node.next_node() == root->child(1)->child(6)->child(1));
      unit_assert (it_node.next_node() == root->child(1)->child(6)->child(2));
      unit_assert (it_node.next_node() == root->child(1)->child(6)->child(3));
      unit_assert (it_node.next_node() == root->child(1)->child(6)->child(4));
      unit_assert (it_node.next_node() == root->child(1)->child(6)->child(5));
      unit_assert (it_node.next_node() == root->child(1)->child(6)->child(6));
      unit_assert (it_node.next_node() == root->child(1)->child(6)->child(7));
      unit_assert (it_node.next_node() == root->child(1)->child(7));
      unit_assert (it_node.next_node() == root->child(2));
      unit_assert (it_node.next_node() == root->child(2)->child(0));
      unit_assert (it_node.next_node() == root->child(2)->child(1));
      unit_assert (it_node.next_node() == root->child(2)->child(2));
      unit_assert (it_node.next_node() == root->child(2)->child(3));
      unit_assert (it_node.next_node() == root->child(2)->child(4));
      unit_assert (it_node.next_node() == root->child(2)->child(5));
      unit_assert (it_node.next_node() == root->child(2)->child(5)->child(0));
      unit_assert (it_node.next_node() == root->child(2)->child(5)->child(1));
      unit_assert (it_node.next_node() == root->child(2)->child(5)->child(2));
      unit_assert (it_node.next_node() == root->child(2)->child(5)->child(3));
      unit_assert (it_node.next_node() == root->child(2)->child(5)->child(4));
      unit_assert (it_node.next_node() == root->child(2)->child(5)->child(5));
      unit_assert (it_node.next_node() == root->child(2)->child(5)->child(6));
      unit_assert (it_node.next_node() == root->child(2)->child(5)->child(7));
      unit_assert (it_node.next_node() == root->child(2)->child(6));
      unit_assert (it_node.next_node() == root->child(2)->child(7));
      unit_assert (it_node.next_node() == root->child(3));
      unit_assert (it_node.next_node() == root->child(3)->child(0));
      unit_assert (it_node.next_node() == root->child(3)->child(1));
      unit_assert (it_node.next_node() == root->child(3)->child(2));
      unit_assert (it_node.next_node() == root->child(3)->child(3));
      unit_assert (it_node.next_node() == root->child(3)->child(4));
      unit_assert (it_node.next_node() == root->child(3)->child(4)->child(0));
      unit_assert (it_node.next_node() == root->child(3)->child(4)->child(1));
      unit_assert (it_node.next_node() == root->child(3)->child(4)->child(2));
      unit_assert (it_node.next_node() == root->child(3)->child(4)->child(3));
      unit_assert (it_node.next_node() == root->child(3)->child(4)->child(4));
      unit_assert (it_node.next_node() == root->child(3)->child(4)->child(5));
      unit_assert (it_node.next_node() == root->child(3)->child(4)->child(6));
      unit_assert (it_node.next_node() == root->child(3)->child(4)->child(7));
      unit_assert (it_node.next_node() == root->child(3)->child(5));
      unit_assert (it_node.next_node() == root->child(3)->child(6));
      unit_assert (it_node.next_node() == root->child(3)->child(7));
      unit_assert (it_node.next_node() == root->child(4));
      unit_assert (it_node.next_node() == root->child(4)->child(0));
      unit_assert (it_node.next_node() == root->child(4)->child(1));
      unit_assert (it_node.next_node() == root->child(4)->child(2));
      unit_assert (it_node.next_node() == root->child(4)->child(3));
      unit_assert (it_node.next_node() == root->child(4)->child(3)->child(0));
      unit_assert (it_node.next_node() == root->child(4)->child(3)->child(1));
      unit_assert (it_node.next_node() == root->child(4)->child(3)->child(2));
      unit_assert (it_node.next_node() == root->child(4)->child(3)->child(3));
      unit_assert (it_node.next_node() == root->child(4)->child(3)->child(4));
      unit_assert (it_node.next_node() == root->child(4)->child(3)->child(5));
      unit_assert (it_node.next_node() == root->child(4)->child(3)->child(6));
      unit_assert (it_node.next_node() == root->child(4)->child(3)->child(7));
      unit_assert (it_node.next_node() == root->child(4)->child(4));
      unit_assert (it_node.next_node() == root->child(4)->child(5));
      unit_assert (it_node.next_node() == root->child(4)->child(6));
      unit_assert (it_node.next_node() == root->child(4)->child(7));
      unit_assert (it_node.next_node() == root->child(5));
      unit_assert (it_node.next_node() == root->child(5)->child(0));
      unit_assert (it_node.next_node() == root->child(5)->child(1));
      unit_assert (it_node.next_node() == root->child(5)->child(2));
      unit_assert (it_node.next_node() == root->child(5)->child(2)->child(0));
      unit_assert (it_node.next_node() == root->child(5)->child(2)->child(1));
      unit_assert (it_node.next_node() == root->child(5)->child(2)->child(2));
      unit_assert (it_node.next_node() == root->child(5)->child(2)->child(3));
      unit_assert (it_node.next_node() == root->child(5)->child(2)->child(4));
      unit_assert (it_node.next_node() == root->child(5)->child(2)->child(5));
      unit_assert (it_node.next_node() == root->child(5)->child(2)->child(6));
      unit_assert (it_node.next_node() == root->child(5)->child(2)->child(7));
      unit_assert (it_node.next_node() == root->child(5)->child(3));
      unit_assert (it_node.next_node() == root->child(5)->child(4));
      unit_assert (it_node.next_node() == root->child(5)->child(5));
      unit_assert (it_node.next_node() == root->child(5)->child(6));
      unit_assert (it_node.next_node() == root->child(5)->child(7));
      unit_assert (it_node.next_node() == root->child(6));
      unit_assert (it_node.next_node() == root->child(6)->child(0));
      unit_assert (it_node.next_node() == root->child(6)->child(1));
      unit_assert (it_node.next_node() == root->child(6)->child(1)->child(0));
      unit_assert (it_node.next_node() == root->child(6)->child(1)->child(1));
      unit_assert (it_node.next_node() == root->child(6)->child(1)->child(2));
      unit_assert (it_node.next_node() == root->child(6)->child(1)->child(3));
      unit_assert (it_node.next_node() == root->child(6)->child(1)->child(4));
      unit_assert (it_node.next_node() == root->child(6)->child(1)->child(5));
      unit_assert (it_node.next_node() == root->child(6)->child(1)->child(6));
      unit_assert (it_node.next_node() == root->child(6)->child(1)->child(7));
      unit_assert (it_node.next_node() == root->child(6)->child(2));
      unit_assert (it_node.next_node() == root->child(6)->child(3));
      unit_assert (it_node.next_node() == root->child(6)->child(4));
      unit_assert (it_node.next_node() == root->child(6)->child(5));
      unit_assert (it_node.next_node() == root->child(6)->child(6));
      unit_assert (it_node.next_node() == root->child(6)->child(7));
      unit_assert (it_node.next_node() == root->child(7));
      unit_assert (it_node.next_node() == root->child(7)->child(0));
      unit_assert (it_node.next_node() == root->child(7)->child(0)->child(0));
      unit_assert (it_node.next_node() == root->child(7)->child(0)->child(1));
      unit_assert (it_node.next_node() == root->child(7)->child(0)->child(2));
      unit_assert (it_node.next_node() == root->child(7)->child(0)->child(3));
      unit_assert (it_node.next_node() == root->child(7)->child(0)->child(4));
      unit_assert (it_node.next_node() == root->child(7)->child(0)->child(5));
      unit_assert (it_node.next_node() == root->child(7)->child(0)->child(6));
      unit_assert (it_node.next_node() == root->child(7)->child(0)->child(7));
      unit_assert (it_node.next_node() == root->child(7)->child(1));
      unit_assert (it_node.next_node() == root->child(7)->child(2));
      unit_assert (it_node.next_node() == root->child(7)->child(3));
      unit_assert (it_node.next_node() == root->child(7)->child(4));
      unit_assert (it_node.next_node() == root->child(7)->child(5));
      unit_assert (it_node.next_node() == root->child(7)->child(6));
      unit_assert (it_node.next_node() == root->child(7)->child(7));
      unit_assert (it_node.next_node() == 0);
      unit_assert (it_node.done()==true);
    }
  }
  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

