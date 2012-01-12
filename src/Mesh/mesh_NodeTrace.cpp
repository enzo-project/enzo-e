// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_NodeTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-01-09
/// @brief    Implementation of the NodeTrace class

#include "mesh.hpp"

//----------------------------------------------------------------------

NodeTrace::NodeTrace(int d, int r) throw ()
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

Node * NodeTrace::trace (Node * node, int ix, int iy, int iz)
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

void NodeTrace::reverse_trace()
{
  unsigned long long ecart[3] = {0,0,0};

  // Reverse bits in trace_[], e.g. trace_[0] = ...00001101 becomes 10110000...

  for (int i=0; i<d_; i++) {

    // Reverse bits byte-by-byte

    char *tr = (char *)&trace_[i];
    char *rt = (char *)&ecart [i];

    *(rt+0) = ((*(tr+7) * 0x80200802ULL) & 0x0884422110ULL) 
      * 0x0101010101ULL >> 32;
    *(rt+1) = ((*(tr+6) * 0x80200802ULL) & 0x0884422110ULL) 
      * 0x0101010101ULL >> 32;
    *(rt+2) = ((*(tr+5) * 0x80200802ULL) & 0x0884422110ULL) 
      * 0x0101010101ULL >> 32;
    *(rt+3) = ((*(tr+4) * 0x80200802ULL) & 0x0884422110ULL)
      * 0x0101010101ULL >> 32;
    *(rt+4) = ((*(tr+3) * 0x80200802ULL) & 0x0884422110ULL) 
      * 0x0101010101ULL >> 32;
    *(rt+5) = ((*(tr+2) * 0x80200802ULL) & 0x0884422110ULL) 
      * 0x0101010101ULL >> 32;
    *(rt+6) = ((*(tr+1) * 0x80200802ULL) & 0x0884422110ULL) 
      * 0x0101010101ULL >> 32;
    *(rt+7) = ((*(tr+0) * 0x80200802ULL) & 0x0884422110ULL)
      * 0x0101010101ULL >> 32;

    ecart[i] = ecart[i] >> 64-level_;


    trace_[i] = ecart[i];
  }

  
}


