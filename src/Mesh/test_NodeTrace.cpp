// See LICENSE_CELLO file for license and copyright information

/// @file     test_NodeTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-02
/// @brief    Test program for the NodeTrace class

#include "main.hpp"
#include "test.hpp"

#include "mesh.hpp"

PARALLEL_MAIN_BEGIN
{

  PARALLEL_INIT;

  unit_init(0,1);

  unit_class("NodeTrace");

  // Create hierarchy

  int d = 2;
  int r = 2;
  int r2d = 4;
  unit_func("Node setup");

  Node *node0, *node1, *node2, *node3, *node4;

  node0 = new Node;
  node0->refine(r2d);
  node1 = node0->child(0);   // position (0,0)
  node1->refine(r2d);
  node2 = node1->child(1);   // position (1,0)
  node2->refine(r2d);
  node3 = node2->child(2);   // position (0,1)
  node3->refine(r2d);
  node4 = node3->child(3);   // position (1,1)

  unit_assert (! node0 -> is_leaf());
  unit_assert (! node1 -> is_leaf());
  unit_assert (! node2 -> is_leaf());
  unit_assert (! node3 -> is_leaf());
  unit_assert (  node4 -> is_leaf());

  //  NodeTrace * NodeTrace = new NodeTrace (d,r);

  unit_func ("NodeTrace");

  NodeTrace * node_trace = new NodeTrace (d,r);
  unit_assert (node_trace != NULL);

  //--------------------------------------------------

  unit_func ("trace()");

  unit_assert(node_trace->trace(node0,0,0) == node1);
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node0,1,0)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node0,0,1)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node0,1,1)->is_leaf());
  node_trace->reset_trace();

  unit_assert(node_trace->trace(node1,0,0)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node1,1,0) == node2);
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node1,0,1)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node1,1,1)->is_leaf());
  node_trace->reset_trace();

  unit_assert(node_trace->trace(node2,0,0)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node2,1,0)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node2,0,1) == node3);
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node2,1,1)->is_leaf());
  node_trace->reset_trace();

  unit_assert(node_trace->trace(node3,0,0)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node3,1,0)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node3,0,1)->is_leaf());
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node3,1,1) == node4); // FAILS
  node_trace->reset_trace();

  unit_assert(node_trace->trace(node4,0,0) == 0);
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node4,1,0) == 0);
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node4,0,1) == 0);
  node_trace->reset_trace();
  unit_assert(node_trace->trace(node4,1,1) == 0);
  node_trace->reset_trace();

  //--------------------------------------------------

  unit_finalize();

  PARALLEL_EXIT;
}

PARALLEL_MAIN_END

