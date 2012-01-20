// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItNode.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-11
/// @brief    Implementation of the ItNode iterator class
//----------------------------------------------------------------------


#include "mesh.hpp"

//----------------------------------------------------------------------

ItNode::ItNode( Tree * tree) throw ()
  : tree_(tree),
    node_trace_(tree->root_node()),
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
    seek_first_leaf_();

  } else {

    seek_next_fork_();

    if (node_trace_.level() == 0) {

      reset_ = true;
      return 0;

    } else {

      seek_next_leaf_();

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
