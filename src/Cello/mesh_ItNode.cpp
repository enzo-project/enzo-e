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
  : lower_level_(0),
    upper_level_(std::numeric_limits<int>::max()),
    tree_(tree),
    node_trace_(tree->root_node()),
    reset_(true)
{
}

//----------------------------------------------------------------------

ItNode::ItNode( Tree * tree, int level) throw ()
  : lower_level_(level),
    upper_level_(level),
    tree_(tree),
    node_trace_(tree->root_node()),
    reset_(true)
{
}

//----------------------------------------------------------------------

ItNode::ItNode( Tree * tree,
		int lower_level,
		int upper_level) throw ()
  : lower_level_(lower_level),
    upper_level_(upper_level),
    tree_(tree),
    node_trace_(tree->root_node()),
    reset_(true)
{
}

//----------------------------------------------------------------------

ItNode::~ItNode() throw ()
{
}

//======================================================================

Node * ItNode::next_leaf() throw()
{
  // if previous was reset, go to first node

  bool in_range, is_leaf;
  do {
    if (reset_) {

      reset_ = false;
      node_trace_.reset();
      seek_next_leaf_();

    } else {

      seek_next_fork_();

      if (node_trace_.level() != 0) {

	seek_next_sibling_();
	seek_next_leaf_();

      } else {

	reset_ = true;
	return 0;

      }
    }
    int level = node_trace_.level();
    in_range = (lower_level_ <= level) && (level <= upper_level_);
    is_leaf  = node_trace_.node()->is_leaf();
  } while ( ! (in_range && is_leaf));

  return node_trace_.node();

}

//----------------------------------------------------------------------

bool ItNode::done () const throw()
{
  return reset_;
}
//----------------------------------------------------------------------
