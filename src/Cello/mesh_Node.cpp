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
  : 
#ifdef CONFIG_USE_CHARM
  have_data_(0),
#endif
  data_(0),
#ifdef CONFIG_USE_CHARM
  size_(0),
#endif
    child_array_(0)
{}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM
void Node::pup (PUP::er &p)
{
  // NOTE: change this function whenever attributes change

  TRACEPUP;
  bool up = p.isUnpacking();

  TRACE1("data_ = %p",data_);
  // int * test = new int[1];
  // *test = 2;
  // int *& test2 = test;
  // TRACE1 ("test2 = %d",*test2);
  // Patch *& data_alias = (Patch *) data_;
  p | have_data_;
  if (have_data_) {
#ifdef REMOVE_PATCH
    if (up) data_ = (void *) new CProxy_CommBlock;
    p | *((CProxy_CommBlock *)data_);
#else /* REMOVE_PATCH */
    if (up) data_ = (void *) new CProxy_Patch;
    p | *((CProxy_Patch *)data_);
#endif
  }
  p | size_;
  if (up) child_array_ = new Node[size_];
  PUParray(p,child_array_,size_);
};
#endif /* CONFIG_USE_CHARM */

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

