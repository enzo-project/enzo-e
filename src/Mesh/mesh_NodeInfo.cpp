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

  if      (r >= 16) bits_ = 4;
  else if (r >= 8)  bits_ = 3;
  else if (r >= 4)  bits_ = 2;
  else if (r >= 2)  bits_ = 1;

  trace_[0] = 0;
  trace_[1] = 0;
  trace_[2] = 0;
}

//----------------------------------------------------------------------

Node * NodeInfo::trace (Node * node, int ix, int iy, int iz)
/// @param node   The node we're tracing from
/// @param ix,iy,iz direction of the child to trace to
/// @return       The child Node we traced to, or NULL if done
{

  int index_child = ix + r_*(iy + r_*iz);

  Node * child = node->child(index_child);

  if (child != NULL) {
    trace_[0] = (trace_[0] << bits_) + ix;
    trace_[1] = (trace_[1] << bits_) + iy;
    trace_[2] = (trace_[2] << bits_) + iz;
  }


  return child;
}

//----------------------------------------------------------------------

void NodeInfo::reverse_trace()
{
  unsigned long long ecart[3] = {0,0,0};
  char *tr[3], *rt[3];
  tr[0]=(char *)&trace_[0];
  rt[0]=(char *)&ecart[0];
  tr[1]=(char *)&trace_[1];
  rt[1]=(char *)&ecart[1];
  tr[2]=(char *)&trace_[2];
  rt[2]=(char *)&ecart[2];

  // Reverse bits in trace_[], e.g. trace_[0] = ...00001101 becomes 10110000...

  *(rt[0]+0) = ((*(tr[0]+7) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[0]+1) = ((*(tr[0]+6) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[0]+2) = ((*(tr[0]+5) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[0]+3) = ((*(tr[0]+4) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[0]+4) = ((*(tr[0]+3) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[0]+5) = ((*(tr[0]+2) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[0]+6) = ((*(tr[0]+1) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[0]+7) = ((*(tr[0]+0) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;

  *(rt[1]+0) = ((*(tr[1]+7) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[1]+1) = ((*(tr[1]+6) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[1]+2) = ((*(tr[1]+5) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[1]+3) = ((*(tr[1]+4) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[1]+4) = ((*(tr[1]+3) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[1]+5) = ((*(tr[1]+2) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[1]+6) = ((*(tr[1]+1) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[1]+7) = ((*(tr[1]+0) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;

  *(rt[2]+0) = ((*(tr[2]+7) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[2]+1) = ((*(tr[2]+6) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[2]+2) = ((*(tr[2]+5) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[2]+3) = ((*(tr[2]+4) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[2]+4) = ((*(tr[2]+3) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[2]+5) = ((*(tr[2]+2) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[2]+6) = ((*(tr[2]+1) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
  *(rt[2]+7) = ((*(tr[2]+0) * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;

  trace_[0] = ecart[0];
  trace_[1] = ecart[1];
  trace_[2] = ecart[2];
}


