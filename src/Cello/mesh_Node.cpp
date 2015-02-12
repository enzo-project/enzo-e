// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Node.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-26

#include "mesh.hpp"

//----------------------------------------------------------------------

Node::Node() throw ()
  : 
  have_data_(0),
  data_(0),
  child_array_()
{}

//----------------------------------------------------------------------

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
    if (up) data_ = (void *) new CProxy_Block;
    p | *((CProxy_Block *)data_);
  }
  p | child_array_;
};

//----------------------------------------------------------------------

Node::~Node() throw ()
{
}

//----------------------------------------------------------------------

int Node::refine (int c)
{
  if (child_array_.size() == 0) {
    child_array_.resize(c);
  } else {
    ERROR ("Node::refine","Cannot refine a Node that has already been refined");
  }
  return c;
}

//----------------------------------------------------------------------

int Node::coarsen (int c)
{
  int count = 0;
  if (child_array_.size() != 0) {
    for (int i=0; i<c; i++) {
      count += child_array_[i].coarsen(c);
    }
    child_array_.resize(0);
    count += c;
  }
  return count;
}

//----------------------------------------------------------------------

Node * Node::child (int k)
{
  return (child_array_.size() != 0) ? &child_array_.at(k) : 0;
}

//======================================================================

