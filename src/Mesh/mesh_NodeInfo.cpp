// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_NodeInfo.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-09
/// @brief    Implementation of the NodeInfo class

#include "mesh.hpp"

//----------------------------------------------------------------------

NodeInfo::NodeInfo(int d, int r) throw ()
  : d_(d), r_(r), level_(0), bits_(0)
{

  if     (16 <= r) bits_ = 4;
  else if (8 <= r) bits_ = 3;
  else if (4 <= r) bits_ = 2;
  else if (2 <= r) bits_ = 1;

  trace_[0] = 0;
  trace_[1] = 0;
  trace_[2] = 0;
}

//----------------------------------------------------------------------

bool NodeInfo::trace (Node * node, int ix, int iy, int iz)
{
  bool done = true;
  // bit-shift left x-bits (1 for octree, 2 for 4^3 tree, etc) and add index

  // if child class, we're done: reverse bits and tag as leaf
  return done;
}

//----------------------------------------------------------------------



