// See LICENSE_CELLO file for license and copyright information

/// @file     test_NodeTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-11
/// @brief    Test program for the NodeTrace class

#include "main.hpp"

#include "test.hpp"
#include "test_mesh.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("NodeTrace");

  Tree * tree = test_tree_22();
  int max_depth = 4;

  unit_func ("NodeTrace");

  NodeTrace * node_trace = new NodeTrace (tree->root_node(),max_depth);

  unit_assert (node_trace != NULL);

  //--------------------------------------------------


  // root node

  unit_func ("level()");
  unit_assert(node_trace->level() == 0);

  unit_func ("node()");
  unit_assert(node_trace->node() == tree->root_node());

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());

  // node 0

  node_trace->push(0);

  unit_func ("level()");
  unit_assert(node_trace->level() == 1);

  unit_func ("node()");
  unit_assert(node_trace->node()  == tree->root_node()->child(0));

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());
  unit_assert(node_trace->node_level(1) == tree->root_node()->child(0));

  unit_func ("index()");
  unit_assert(node_trace->index() == 0);

  unit_func ("index_level()");
  unit_assert(node_trace->index_level(1) == 0);

  node_trace->reset();

  // node 1

  node_trace->push(1);

  unit_func ("level()");
  unit_assert(node_trace->level() == 1);

  unit_func ("node()");
  unit_assert(node_trace->node()  == tree->root_node()->child(1));

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());
  unit_assert(node_trace->node_level(1) == tree->root_node()->child(1));

  unit_func ("index()");
  unit_assert(node_trace->index() == 1);

  unit_func ("index_level()");
  unit_assert(node_trace->index_level(1) == 1);

  node_trace->reset();

  // node 3

  node_trace->push(3);

  unit_func ("level()");
  unit_assert(node_trace->level() == 1);

  unit_func ("node()");
  unit_assert(node_trace->node()  == tree->root_node()->child(3));

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());
  unit_assert(node_trace->node_level(1) == tree->root_node()->child(3));

  unit_func ("index()");
  unit_assert(node_trace->index() == 3);

  unit_func ("index_level()");
  unit_assert(node_trace->index_level(1) == 3);

  node_trace->reset();

  // node 01

  node_trace->push(0);
  node_trace->push(1);

  unit_func ("level()");
  unit_assert(node_trace->level() == 2);

  unit_func ("node()");
  unit_assert(node_trace->node()  == tree->root_node()->child(0)->child(1));

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());
  unit_assert(node_trace->node_level(1) == tree->root_node()->child(0));
  unit_assert
    (node_trace->node_level(2) == tree->root_node()->child(0)->child(1));

  unit_func ("index()");
  unit_assert(node_trace->index() == 1);

  unit_func ("index_level()");
  unit_assert(node_trace->index_level(1) == 0);
  unit_assert(node_trace->index_level(2) == 1);

  node_trace->reset();

  // node 02

  node_trace->push(0);
  node_trace->push(2);

  unit_func ("level()");
  unit_assert(node_trace->level() == 2);

  unit_func ("node()");
  unit_assert(node_trace->node()  == tree->root_node()->child(0)->child(2));

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());
  unit_assert(node_trace->node_level(1) == tree->root_node()->child(0));
  unit_assert
    (node_trace->node_level(2) == tree->root_node()->child(0)->child(2));

  unit_func ("index()");
  unit_assert(node_trace->index() == 2);

  unit_func ("index_level()");
  unit_assert(node_trace->index_level(1) == 0);
  unit_assert(node_trace->index_level(2) == 2);

  node_trace->reset();


  // node 30

  node_trace->push(3);
  node_trace->push(0);

  unit_func ("level()");
  unit_assert(node_trace->level() == 2);

  unit_func ("node()");
  unit_assert(node_trace->node()  == tree->root_node()->child(3)->child(0));

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());
  unit_assert(node_trace->node_level(1) == tree->root_node()->child(3));
  unit_assert
    (node_trace->node_level(2) == tree->root_node()->child(3)->child(0));

  unit_func ("index()");
  unit_assert(node_trace->index() == 0);

  unit_func ("index_level()");
  unit_assert(node_trace->index_level(1) == 3);
  unit_assert(node_trace->index_level(2) == 0);

  node_trace->reset();


  // node 33

  node_trace->push(3);
  node_trace->push(3);

  unit_func ("level()");
  unit_assert(node_trace->level() == 2);

  unit_func ("node()");
  unit_assert(node_trace->node()  == tree->root_node()->child(3)->child(3));

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());
  unit_assert(node_trace->node_level(1) == tree->root_node()->child(3));
  unit_assert
    (node_trace->node_level(2) == tree->root_node()->child(3)->child(3));

  unit_func ("index()");
  unit_assert(node_trace->index() == 3);

  unit_func ("index_level()");
  unit_assert(node_trace->index_level(1) == 3);
  unit_assert(node_trace->index_level(2) == 3);

  node_trace->reset();


  // node 012

  node_trace->push(0);
  node_trace->push(1);
  node_trace->push(2);

  unit_func ("level()");
  unit_assert(node_trace->level() == 3);

  unit_func ("node()");
  unit_assert(node_trace->node() == 
	      tree->root_node()->child(0)->child(1)->child(2));

  unit_func ("node_level()");
  unit_assert(node_trace->node_level(0) == tree->root_node());
  unit_assert
    (node_trace->node_level(1) == tree->root_node()->child(0));
  unit_assert
    (node_trace->node_level(2) == tree->root_node()->child(0)->child(1));
  unit_assert
    (node_trace->node_level(3) == tree->root_node()->child(0)->child(1)->child(2));

  unit_func ("index()");
  unit_assert(node_trace->index() == 2);

  unit_func ("index_level()");
  unit_assert(node_trace->index_level(1) == 0);
  unit_assert(node_trace->index_level(2) == 1);
  unit_assert(node_trace->index_level(3) == 2);

  node_trace->reset();


  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

