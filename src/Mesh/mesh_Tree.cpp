// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Tree.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "mesh.hpp"

//----------------------------------------------------------------------

Tree::Tree(int d, int r) throw ()
  : d_(d), r_(r), 
    root_(new Node), 
    num_nodes_(1)
{
}

//----------------------------------------------------------------------

Tree::~Tree() throw ()
{
  delete root_;
}

//----------------------------------------------------------------------

void Tree::refine_node (NodeTrace * node_trace)
{
  int    level = node_trace->level();
  Node * node  = node_trace->node();

  int r2d = 1;
  if (d_>=1) r2d *= r_;
  if (d_>=2) r2d *= r_;
  if (d_>=3) r2d *= r_;
  node->refine(r2d);
  
  num_nodes_ += r2d;
}

//----------------------------------------------------------------------

void Tree::delete_node (NodeTrace * node_trace)
{
}

//----------------------------------------------------------------------

void Tree::balance_node (NodeTrace * node_trace)
{
}

//----------------------------------------------------------------------

NodeTrace * Tree::node_neighbor (NodeTrace * node_trace) const
{
  return NULL;
}

//----------------------------------------------------------------------

NodeTrace * Tree::node_parent (NodeTrace * node_trace) const
{
  return NULL;
}

//----------------------------------------------------------------------

NodeTrace * Tree::node_child (NodeTrace * node_trace, int k) const
{
  return NULL;
}

//======================================================================

