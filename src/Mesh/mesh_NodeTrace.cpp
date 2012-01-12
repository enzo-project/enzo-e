// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_NodeTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-09
/// @brief    Implementation of the NodeTrace class

#include "mesh.hpp"

//----------------------------------------------------------------------

NodeTrace::NodeTrace(Node * root, int max_levels) throw ()
  : max_levels_(max_levels),
    level_(0)
{
  TRACE1("max_levels = %d",max_levels_);
  index_ = new int  [ max_levels ];
  index_[0] = -1;
  node_  = new Node * [ max_levels ];
  node_[0] = root;
}

//----------------------------------------------------------------------

NodeTrace::~NodeTrace() throw ()
{
  delete index_;
  delete node_;
}

//----------------------------------------------------------------------

Node * NodeTrace::push (int index)
/// @param node   The node we're tracing from
/// @param index  Index of the node's child we're tracing to
/// @return       The child Node we traced to, or NULL if node is a leaf
{

  ASSERT1 ("NodeTrace::push",
	   "Trying to push more nodes than the %d allowed",
	   max_levels_,
	   level_ < max_levels_ - 1);

  Node * child = node_[level_]->child(index);

  if (child != NULL) {
    ++level_;
    index_[level_] = index;
    node_[level_]  = child;
  }

  return child;
}

//----------------------------------------------------------------------

Node * NodeTrace::pop ()
{
  ASSERT ("NodeTrace::pop", 
	  "Trying to pop with level_ == 0", 
	  level_ > 0);
  --level_;
  return node_[level_];
}

//----------------------------------------------------------------------

Node * NodeTrace::node() const
{ 
  return node_[level_];
};

//----------------------------------------------------------------------

Node * NodeTrace::node_level(int level) const
{ 
  ASSERT2 ("NodeTrace::node", 
	   "input level = %d is not between 0 and node level %d",
	   level,level_,
	   0 <= level && level <= level_);
  return node_[level];
};

//----------------------------------------------------------------------

int NodeTrace::index() const
{
  return index_[level_];
}

//----------------------------------------------------------------------

int NodeTrace::index_level(int level) const
{
  ASSERT2 ("NodeTrace::index",
	   "input level = %d is not between 0 and node level %d",
	   level,level_,
	   0 <= level && level <= level_);
  return index_[level];
}
