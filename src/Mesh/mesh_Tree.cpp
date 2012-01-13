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
  c_ = 1;
  if (d_>=1) c_ *= r_;
  if (d_>=2) c_ *= r_;
  if (d_>=3) c_ *= r_;
}

//----------------------------------------------------------------------

Tree::~Tree() throw ()
{
  delete root_;
}

//----------------------------------------------------------------------

void Tree::refine_node (NodeTrace * node_trace)
{
  node_trace->node()->refine(c_);
  
  num_nodes_ += c_;
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

