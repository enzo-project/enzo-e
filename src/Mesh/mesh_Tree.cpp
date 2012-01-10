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
    num_nodes_(1), 
    max_level_(0)
{
}

//----------------------------------------------------------------------

Tree::~Tree() throw ()
{
  delete root_;
}

//----------------------------------------------------------------------

void Tree::refine_node (Node * node)
{
  int r2d = 1;
  if (d_>=1) r2d *= r_;
  if (d_>=2) r2d *= r_;
  if (d_>=3) r2d *= r_;
  node->refine(r2d);
  num_nodes_ += r2d;
}

//----------------------------------------------------------------------

void Tree::delete_node (Node * node)
{
}

//----------------------------------------------------------------------

void Tree::balance_node (Node * node)
{
}

//----------------------------------------------------------------------

Node * Tree::node_neighbor (Node * node) const
{
  return NULL;
}

//----------------------------------------------------------------------

Node * Tree::node_parent (Node * node) const
{
  return NULL;
}

//----------------------------------------------------------------------

Node * Tree::node_child (Node * node, int k) const
{
  return NULL;
}

//======================================================================

