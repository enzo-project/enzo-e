// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Node.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    
///
/// 

#include "mesh.hpp"

//----------------------------------------------------------------------

Node::Node() throw ()
  : data_(0),child_(0)
{}

//----------------------------------------------------------------------

Node::~Node() throw ()
{
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
void Node::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change
  p | *((Patch * )data_);
  p | *child_;
}
#endif

//----------------------------------------------------------------------

int Node::refine (int c)
{
  if (child_ == 0) {
    child_ = new Node [c];
  } else {
    ERROR ("Node::refine","Cannot refine a Node that has already been refined");
  }
  return c;
}

//----------------------------------------------------------------------

int Node::coarsen (int c)
{
  int count = 0;
  if (child_ != 0) {
    for (int i=0; i<c; i++) {
      count += child_[i].coarsen(c);
    }
    delete [] child_;
    child_ = 0;
    count += c;
  }
  return count;
}

//----------------------------------------------------------------------

Node * Node::child (int k) const
{
  return (child_ != 0) ? &child_[k] : 0;
}

//======================================================================

