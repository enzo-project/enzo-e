// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-19
/// @todo     Clean up interface
/// @brief    Implementation of the Layout class

//----------------------------------------------------------------------

#include "cello.hpp"

#include "parallel.hpp"

//----------------------------------------------------------------------

Layout::Layout(int nbx, int nby, int nbz) throw()
  : process_first_ (0),
    process_count_ (1)
{
  block_count_[0] = nbx;
  block_count_[1] = nby;
  block_count_[2] = nbz;

}

//----------------------------------------------------------------------

void Layout::set_process_range(int process_first, int process_count) throw()
{
  process_first_ = process_first;
  process_count_  = process_count;
}

//----------------------------------------------------------------------

int Layout::block_count (int *nbx, int *nby, int *nbz) throw()
{ 
  *nbx = block_count_[0];
  *nby = block_count_[1];
  *nbz = block_count_[2];
  return (*nbx)*(*nby)*(*nbz);
}

//----------------------------------------------------------------------

void Layout::process_range(int * process_first, int * process_count) throw()
{
  *process_first = process_first_;
  *process_count = process_count_;
}

//----------------------------------------------------------------------

int Layout::local_count (int ip) throw()
{
  int ip0 = ip - process_first_;

  if (0 <= ip0 && ip0 < process_count_) {

    int block_count = block_count_[0] * block_count_[1] * block_count_[2];

    return (ip0+1)*block_count/process_count_ 
      -     ip0   *block_count/process_count_;

  } else {
    return 0;
  }
}

//----------------------------------------------------------------------

bool Layout::is_local (int ip, int ibx, int iby, int ibz) throw()
{
  WARNING("Layout::is_local",
	  "index_first_local is never != for current test code");
  int block_count = block_count_[0] * block_count_[1] * block_count_[2];
  int ib = block_index(ibx,iby,ibz);
  int index_first_local = (ip-process_first_)*block_count/process_count_;
  return 0 <= ib-index_first_local && ib-index_first_local < local_count(ip);
}

//----------------------------------------------------------------------

int Layout::global_index (int ip, int ib) throw()
{
  int block_count = block_count_[0] * block_count_[1] * block_count_[2];
  return ib + (ip-process_first_)*block_count/process_count_;
}

//----------------------------------------------------------------------

int Layout::process (int ib)  throw()
{
  int block_count = block_count_[0] * block_count_[1] * block_count_[2];
  if (0 <= ib && ib < block_count) {
    int ip0 = process_count_*ib / block_count;
      return process_first_ + ip0;
  } else {
    return PROCESS_NULL;
  }
}

//----------------------------------------------------------------------

int Layout::process (int ibx, int iby, int ibz)  throw()
{
  int nb = block_index (ibx,iby,ibz);
  return process(nb);
}
//----------------------------------------------------------------------

int Layout::block_index (int ibx, int iby, int ibz) throw()
{
  return ibx + block_count_[0]*(iby + block_count_[1]*ibz);
}

//----------------------------------------------------------------------

void Layout::block_indices (int ib, int * ibx, int * iby, int * ibz) throw()
{
  *ibx = ib % block_count_[0];
  *iby = (ib / block_count_[0]) % block_count_[1];
  *ibz = ib / (block_count_[0]*block_count_[1]);
}
