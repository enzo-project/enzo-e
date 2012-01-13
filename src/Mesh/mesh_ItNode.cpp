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
    node_trace_(tree->root_node(),max_depth),
    reset_(true)
{
}

//----------------------------------------------------------------------

ItNode::~ItNode() throw ()
{
}

//======================================================================

Node * ItNode::operator++ () throw()
{
  // if previous was reset, go to first node

  if (reset_) {
    reset_ = false;
    node_trace_.reset();
    while ( ! node_trace_.node()->is_leaf() ) {
      node_trace_.push(0);
    }

  } else {

    while ((node_trace_.level() > 0) && 
	   (node_trace_.index() + 1 == tree_->num_children()) ) {
      node_trace_.pop();
    }

    if (node_trace_.level() == 0) {

      reset_ = true;
      return 0;

    } else {

      int index = node_trace_.index();
      node_trace_.pop();
      node_trace_.push(index + 1);
      while ( ! node_trace_.node()->is_leaf() ) {
	node_trace_.push(0);
      }
    }
  }

  return node_trace_.node();

}
//----------------------------------------------------------------------

bool ItNode::done () const throw()
{
  return reset_;
}
//----------------------------------------------------------------------
