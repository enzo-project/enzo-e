// See LICENSE_CELLO file for license and copyright information

//----------------------------------------------------------------------
/// @file     mesh_ItNode.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-11
/// @brief    Implementation of the ItNode iterator class
//----------------------------------------------------------------------


#include "mesh.hpp"

//----------------------------------------------------------------------

ItNode::ItNode( Tree * tree ) throw ()
  : tree_(tree),
    node_(0),
    node_trace_(new NodeTrace(tree->dimension(),tree->refinement()))
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
  }
}
//----------------------------------------------------------------------

bool ItNode::done () const throw()
{
}
//----------------------------------------------------------------------
