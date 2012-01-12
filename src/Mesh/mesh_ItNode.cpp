// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItNode.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-11
/// @brief    Implementation of the ItNode iterator class
//----------------------------------------------------------------------


#include "mesh.hpp"

//----------------------------------------------------------------------

ItNode::ItNode( Tree * tree , int max_depth) throw ()
  : tree_(tree),
    node_(0),
    node_trace_(new NodeTrace(tree->root_node(),max_depth))
{
}

//----------------------------------------------------------------------

ItNode::~ItNode() throw ()
{
  delete node_trace_;
}

//======================================================================

Node * ItNode::operator++ () throw()
{
  if (node_ == 0) {
    node_ = tree_->root_node();
    node_trace_->reset();
  } else {
    // if node is not a leaf, 
    if (! node_->is_leaf()) {
      
      
    } else {
    }

  }
  return node_;
}
//----------------------------------------------------------------------

bool ItNode::done () const throw()
{
}
//----------------------------------------------------------------------
