// $Id$
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the Layout class

//----------------------------------------------------------------------

#include "error.hpp"
#include "parallel.hpp"

//----------------------------------------------------------------------

Layout::Layout() throw()
  : process_first_ (0),
    process_count_ (1)
{

  for (int i=0; i<3; i++) {
    block_count_[i] = 1;
  }
}

//----------------------------------------------------------------------

void Layout::set_range(int process_first, int process_count) throw()
{
  process_first_ = process_first;
  process_count_  = process_count;
}

//----------------------------------------------------------------------

void Layout::set_blocks(int nb0, int nb1, int nb2) throw()
{
  block_count_[0] = nb0;
  block_count_[1] = nb1;
  block_count_[2] = nb2; 
}

//----------------------------------------------------------------------

int Layout::blocks (int *nb0, int *nb1, int *nb2) throw()
{ 
  *nb0 = block_count_[0];
  *nb1 = block_count_[1];
  *nb2 = block_count_[2];
  return (*nb0)*(*nb1)*(*nb2);
}

//----------------------------------------------------------------------

void Layout::range(int * process_first, int * process_count) throw()
{
  *process_first = process_first_;
  *process_count  = process_count_;
}

//----------------------------------------------------------------------

int Layout::process (int ibx, int iby, int ibz)  throw()
{
  int ib = ibx + block_count_[0]*(iby + block_count_[1]*ibz);
  int nb = block_count_[0] * block_count_[1] * block_count_[2];

  return process_first_ + process_count_*ib / nb;
    
}
