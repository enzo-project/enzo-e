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
  index_ = new int  [ max_levels ];
  index_[0] = -1;
  node_  = new Node * [ max_levels ];
  node_[0] = root;
}

//----------------------------------------------------------------------

NodeTrace::NodeTrace(const NodeTrace & node_trace) throw()
{
  copy_(node_trace);
}

//----------------------------------------------------------------------

NodeTrace & NodeTrace::operator=(const NodeTrace &node_trace) throw()
{
  if (this != &node_trace) {
    copy_(node_trace);
  }
  return *this;
}

//----------------------------------------------------------------------

void NodeTrace::copy_(const NodeTrace & node_trace)
{
  max_levels_ = node_trace.max_levels_;
  level_      = node_trace.level_;
  index_ = new int    [ max_levels_ ];
  for (int i=0; i<max_levels_; i++) {
    index_[i] = node_trace.index_[i];
  }
  node_  = new Node * [ max_levels_ ];
  for (int i=0; i<max_levels_; i++) {
    node_[i] = node_trace.node_[i];
  }
}

//----------------------------------------------------------------------

NodeTrace::~NodeTrace() throw ()
{
  delete index_;
  delete node_;
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
