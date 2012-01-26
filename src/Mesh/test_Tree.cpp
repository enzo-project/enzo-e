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
  unit_func("node_child ()");

  NodeTrace node_trace (tree->root_node());

  NodeTrace child_trace (tree->node_child(node_trace,0));
  unit_assert(*((int *)child_trace.node()->data()) == 1);

  child_trace = tree->node_child(node_trace,1);
  unit_assert(*((int *)child_trace.node()->data()) == 2);

  child_trace = tree->node_child(node_trace,2);
  unit_assert(*((int *)child_trace.node()->data()) == 3);

  child_trace = tree->node_child(node_trace,3);
  unit_assert(*((int *)child_trace.node()->data()) == 4);

  child_trace =  tree->node_child(tree->node_child(node_trace,0),0);
  unit_assert(*((int *)child_trace.node()->data()) == 5);

  child_trace =  tree->node_child(tree->node_child(node_trace,0),1);
  unit_assert(*((int *)child_trace.node()->data()) == 6);

  child_trace =  tree->node_child(tree->node_child(node_trace,0),2);
  unit_assert(*((int *)child_trace.node()->data()) == 7);

  child_trace =  tree->node_child(tree->node_child(node_trace,0),3);
  unit_assert(*((int *)child_trace.node()->data()) == 8);

  child_trace =  
    tree->node_child 
    ( tree->node_child 
      ( tree->node_child ( node_trace, 3), 0), 2);

  unit_assert(*((int *)child_trace.node()->data()) == 27);

  //--------------------------------------------------

  unit_func("node_parent ()");

  NodeTrace parent_trace (node_trace);

  child_trace = tree->node_child(node_trace,0);
  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);

  child_trace = tree->node_child(node_trace,1);
  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);

  child_trace = tree->node_child(node_trace,2);
  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);

  child_trace = tree->node_child(node_trace,3);
  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);

  child_trace =  tree->node_child(tree->node_child(node_trace,0),0);
  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 1);
  parent_trace = tree->node_parent(parent_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);

  child_trace =  tree->node_child(tree->node_child(node_trace,0),1);
  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 1);
  parent_trace = tree->node_parent(parent_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);

  child_trace =  tree->node_child(tree->node_child(node_trace,0),2);
  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 1);
  parent_trace = tree->node_parent(parent_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);

  child_trace =  tree->node_child(tree->node_child(node_trace,0),3);
  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 1);
  parent_trace = tree->node_parent(parent_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);

  child_trace =  
    tree->node_child 
    ( tree->node_child 
      ( tree->node_child ( node_trace, 3), 0), 2);

  parent_trace = tree->node_parent(child_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 21);
  parent_trace = tree->node_parent(parent_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 4);
  parent_trace = tree->node_parent(parent_trace);
  unit_assert(*((int *)parent_trace.node()->data()) == 0);
  


  //--------------------------------------------------
  //--------------------------------------------------
  unit_func("refine_node () (tested in test_ItNode)");
  unit_assert (true);

  //--------------------------------------------------
  unit_func("coarsen_node ()");

  {

    unit_assert(tree->root_node()->child(1)->is_leaf() == false);

    NodeTrace node_trace (tree->root_node());
  
    node_trace.push(1);

    tree->coarsen_node(node_trace);

    unit_assert(tree->num_nodes()==29);
    unit_assert(tree->root_node()->child(1)->is_leaf() == true);

    unit_assert(tree->root_node()->child(0)->is_leaf() == false);

    node_trace.pop();
    node_trace.push(0);
    tree->coarsen_node(node_trace);

    unit_assert(tree->num_nodes()==17);
    unit_assert(tree->root_node()->child(0)->is_leaf() == true);
  
    delete tree;
  }

  //--------------------------------------------------
  unit_func("node_neighbor ()");

  {
    tree = test_tree_22();

    NodeTrace node_trace (tree->root_node());
    NodeTrace trace16(tree->root_node()), neighbor16(tree->root_node());
    trace16.push(0);
    trace16.push(2);
    trace16.push(3);
    tree->node_neighbor(trace16,&neighbor16,0,+1);

    int * data;

    NodeTrace neighbor(node_trace);

    // left neighbor of node 9 is node 5
    node_trace.push(0);
    node_trace.push(1);
    node_trace.push(0);

    unit_assert(tree->node_neighbor(node_trace,&neighbor,-1,0,0) == true);
    data = (int *)neighbor.node()->data();
    unit_assert (*data == 5);

    // right neighbor of node 5 is node 6 (which has children at a finer level)

    node_trace.reset();

    node_trace.push(0);
    node_trace.push(0);

    unit_assert(tree->node_neighbor(node_trace,&neighbor,+1,0,0) == true);
    data = (int *)neighbor.node()->data();
    unit_assert (*data == 6);

    // upper neighbor of node 17 is node 19

    node_trace.reset();

    node_trace.push(1);
    node_trace.push(0);

    data = (int *)node_trace.node()->data();
    unit_assert(tree->node_neighbor(node_trace,&neighbor,0,+1,0) == true);
    data = (int *)neighbor.node()->data();
    unit_assert (*data == 19);

    // lower neighbor of node 26 is node 19

    node_trace.reset();

    node_trace.push(3);
    node_trace.push(0);
    node_trace.push(1);

    data = (int *)node_trace.node()->data();
    unit_assert(tree->node_neighbor(node_trace,&neighbor,0,-1,0) == true);
    data = (int *)neighbor.node()->data();
    unit_assert (*data == 19);
  }


  //--------------------------------------------------
    unit_func("balance_node ()");

  {
    int d = 2;
    int r = 2;
    tree = new Tree (d,r);

    int max_depth=12;
    int width  = 1024;
    int height = 1024;

    int nx,ny;

    int * levels = create_levels_from_image 
      ("input/test_balance.png",&nx,&ny,max_depth);

    create_tree_from_levels (tree, levels, nx, ny);

    create_image_from_tree (tree, "test_tree_initial.png",width,height);

    int count_level[max_depth];
    for (int level = 0; level < max_depth; level ++)
    {
      count_level[level] = 0;
      ItNode it_node(tree,level);
      while (it_node.next_leaf()) {
	++count_level[level];
      }
    }

    int count_tree[max_depth];
    for (int level = 0; level < max_depth; level ++)
    {
      count_tree[level] = 0;
      ItNode it_node(tree,0,level);
      while (it_node.next_leaf()) {
	++count_tree[level];
      }
    }

    unit_func("ItNode levels");
    for (int level=1; level < max_depth; level++) {
      unit_assert(count_tree[level] - count_tree[level-1] ==
		  count_level[level]);
    }
    unit_assert(count_tree[max_depth-1] - 1 == 
		(tree->num_nodes() - 1) * 3 / 4);

    printf ("tree initial num_nodes = %d\n",tree->num_nodes());
    printf ("tree initial max_level = %d\n",tree->max_level());


    //--------------------------------------------------
    // Balance tree
    //--------------------------------------------------
    Timer timer;

    timer.start();

    tree->balance();

    printf ("time to balance = %f s\n",timer.value());

    printf ("tree balanced num_nodes = %d\n",tree->num_nodes());
    printf ("tree balanced max_level = %d\n",tree->max_level());

    create_image_from_tree (tree, "test_tree_balanced.png",width,height);
    
    //--------------------------------------------------
    // Coalesce tree nodes
    //--------------------------------------------------

    timer.clear();
    timer.start();
    tree->coalesce();

    printf ("time to coalesce = %f s\n",timer.value());

    printf ("tree coalesced num_nodes = %d\n",tree->num_nodes());
    printf ("tree coalesced max_level = %d\n",tree->max_level());

    create_image_from_tree (tree, "test_tree_coalesced.png",width,height);

    delete tree;
  }

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

