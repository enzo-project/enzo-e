// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_test.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-12
/// @brief    Functions for Mesh test programs

#include "mesh.hpp"

//----------------------------------------------------------------------

Tree * test_tree_22()
{

  // +-------+---+-+-+
  // |       |   + + +
  // |       +-+-+-+-+
  // |       + + +   |
  // +-+-+---+-+-+---+
  // + + +   |   |   |
  // +-+-+-+-+---+---|
  // |   + + +   |   |
  // +---+-+-+-------+
  //
  // Refined nodes: root, 0 1 3 01 02 30 33

  int d = 2;
  int r = 2;

  Tree * tree = new Tree (d,r);
  int max_levels = 4;

  NodeTrace * node_trace = new NodeTrace (tree->root_node(),max_levels);
  
  // root
  tree->refine_node (node_trace);
  // 0
  node_trace->push(0);
  tree->refine_node (node_trace);
  // 01
  node_trace->push(1);
  tree->refine_node (node_trace);
  // 02
  node_trace->pop();
  node_trace->push(2);
  tree->refine_node (node_trace);
  // 1
  node_trace->pop();
  node_trace->pop();
  node_trace->push(1);
  tree->refine_node (node_trace);
  // 3
  node_trace->pop();
  node_trace->push(3);
  tree->refine_node (node_trace);
  // 30
  node_trace->push(0);
  tree->refine_node (node_trace);
  // 33
  node_trace->pop();
  node_trace->push(3);
  tree->refine_node (node_trace);

  delete node_trace;

  return tree;
}

Tree * test_tree_32()
{

  // 3D version of the following;
  // +-------+---+-+-+
  // |   |   |   |   |
  // +---+-+-+-+-+---+
  // |   +-+-+-+-+   |
  // +---+-+-+-+-+---+
  // |   +-+-+-+-+   |
  // +---+-+-+-+-+---+
  // |   |   |   |   |
  // +---+-+-+-------+
  //
  // Refined nodes: root, 0 07 1 16 2 25 3 34 4 43 5 52 6 61 7 70

  int d = 3;
  int r = 2;

  Tree * tree = new Tree (d,r);
  int max_levels = 4;

  NodeTrace * node_trace = new NodeTrace (tree->root_node(),max_levels);

  // root
  tree->refine_node (node_trace);

  // 0
  node_trace->push(0);
  tree->refine_node (node_trace);
  // 07
  node_trace->push(7);
  tree->refine_node (node_trace);
  // 1
  node_trace->pop();
  node_trace->pop();
  node_trace->push(1);
  tree->refine_node (node_trace);
  // 16
  node_trace->push(6);
  tree->refine_node (node_trace);
  // 2
  node_trace->pop();
  node_trace->pop();
  node_trace->push(2);
  tree->refine_node (node_trace);
  // 25
  node_trace->push(5);
  tree->refine_node (node_trace);
  // 3
  node_trace->pop();
  node_trace->pop();
  node_trace->push(3);
  tree->refine_node (node_trace);
  // 34
  node_trace->push(4);
  tree->refine_node (node_trace);
  // 4
  node_trace->pop();
  node_trace->pop();
  node_trace->push(4);
  tree->refine_node (node_trace);
  // 43
  node_trace->push(3);
  tree->refine_node (node_trace);
  // 5
  node_trace->pop();
  node_trace->pop();
  node_trace->push(5);
  tree->refine_node (node_trace);
  // 52
  node_trace->push(2);
  tree->refine_node (node_trace);
  // 6
  node_trace->pop();
  node_trace->pop();
  node_trace->push(6);
  tree->refine_node (node_trace);
  // 61
  node_trace->push(1);
  tree->refine_node (node_trace);
  // 7
  node_trace->pop();
  node_trace->pop();
  node_trace->push(7);
  tree->refine_node (node_trace);
  // 70
  node_trace->push(0);
  tree->refine_node (node_trace);

  delete node_trace;

  return tree;
}
