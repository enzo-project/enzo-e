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
  : data_(0),
#ifdef CONFIG_USE_CHARM
    size_(0),
#endif
    child_array_(0)
{}

//----------------------------------------------------------------------

Node::~Node() throw ()
{
}

//----------------------------------------------------------------------

int Node::refine (int c)
{
  if (child_array_ == 0) {
#ifdef CONFIG_USE_CHARM
    size_ = c;
#endif
    child_array_ = new Node [c];
  } else {
    ERROR ("Node::refine","Cannot refine a Node that has already been refined");
  }
  return c;
}

//----------------------------------------------------------------------

int Node::coarsen (int c)
{
  int count = 0;
  if (child_array_ != 0) {
    for (int i=0; i<c; i++) {
      count += child_array_[i].coarsen(c);
    }
    delete [] child_array_;
#ifdef CONFIG_USE_CHARM
    size_ = 0;
#endif
    child_array_ = 0;
    count += c;
  }
  return count;
}

//----------------------------------------------------------------------

Node * Node::child (int k) const
{
  return (child_array_ != 0) ? &child_array_[k] : 0;
}

//======================================================================

