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

  delete tree;

  tree = new Tree (2,2);

  unit_func("balance_node ()");

  int nx,ny;
  int * levels = create_levels_from_image ("test_balance.png",&nx,&ny,max_depth);

  int counter=0;
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      int i = ix + nx*iy;
      double x = 1.0*ix / nx;
      double y = 1.0*iy / ny;
      NodeTrace node_trace (tree->root_node(), max_depth);
      int a = levels[i];
      while (--a > 0) {
	int rx = x > 0.5;
	int ry = y > 0.5;
	int r = rx + 2*ry;
	Node * node = node_trace.node();
	if (node->is_leaf()) {
	  tree->refine_node (&node_trace);
	}
	node_trace.push(r);

	int * d = new int;
	*d = counter++;
	node_trace.node()->set_data (d);

	rx *= 2.0;
	ry *= 2.0;
	if (rx > 1.0) rx -= 1.0;
	if (ry > 1.0) ry -= 1.0;
      }
    }
  }

  create_image_from_tree (tree, "test_balance_out.png",512,512,max_depth);

  printf ("counter = %d\n",counter);
  printf ("tree num_nodes = %d\n",tree->num_nodes());

  ItNode it_node(tree, max_depth);
  while (Node * node = ++it_node) {
    int * d = (int *)node->data();
    if (d) printf ("d=%d\n",*d);
  }

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

