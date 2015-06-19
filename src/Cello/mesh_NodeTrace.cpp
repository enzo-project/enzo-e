// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_NodeTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-09
/// @brief    Implementation of the NodeTrace class

#include "mesh.hpp"

//----------------------------------------------------------------------

NodeTrace::NodeTrace(Node * root) throw ()
  : level_(0)
{
  index_.push_back(-1);
  node_.push_back(root);
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
  level_  = node_trace.level_;
  index_  = node_trace.index_;
  node_  = node_trace.node_;
}

//----------------------------------------------------------------------

NodeTrace::~NodeTrace() throw ()
{
}

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
