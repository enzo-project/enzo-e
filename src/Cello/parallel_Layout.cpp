// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_Layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2010-04-19
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

int Layout::block_count (int *nbx, int *nby, int *nbz) const throw()
{ 
  *nbx = block_count_[0];
  *nby = block_count_[1];
  *nbz = block_count_[2];
  return (*nbx)*(*nby)*(*nbz);
}

//----------------------------------------------------------------------

void Layout::process_range
(int * process_first, int * process_count) const throw()
{
  if (process_first) (*process_first) = process_first_;
  if (process_count) (*process_count) = process_count_;
}

//----------------------------------------------------------------------

int Layout::local_count (int ip) const throw()
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

bool Layout::is_local (int ip, int ibx, int iby, int ibz) const throw()
{
  WARNING("Layout::is_local",
	  "index_first_local is never != for current test code");
  int block_count = block_count_[0] * block_count_[1] * block_count_[2];
  int ib = block_index(ibx,iby,ibz);
  int index_first_local = (ip-process_first_)*block_count/process_count_;
  return 0 <= ib-index_first_local && ib-index_first_local < local_count(ip);
}

//----------------------------------------------------------------------

int Layout::global_index (int ip, int ib) const throw()
{
  int block_count = block_count_[0] * block_count_[1] * block_count_[2];
  return ib + (ip-process_first_)*block_count/process_count_;
}

//----------------------------------------------------------------------

int Layout::process (int ib)  const throw()
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

int Layout::process3 (int ibx, int iby, int ibz)  const throw()
{
  // Force block indices to be in range for periodic b.c.
  while (ibx >= block_count_[0]) ibx -= block_count_[0];
  while (ibx < 0)                ibx += block_count_[0];
  while (iby >= block_count_[1]) iby -= block_count_[1];
  while (iby < 0)                iby += block_count_[1];
  while (ibz >= block_count_[2]) ibz -= block_count_[2];
  while (ibz < 0)                ibz += block_count_[2];
  int nb = block_index (ibx,iby,ibz);
  return process(nb);
}
//----------------------------------------------------------------------

int Layout::block_index (int ibx, int iby, int ibz) const throw()
{
  return ibx + block_count_[0]*(iby + block_count_[1]*ibz);
}

//----------------------------------------------------------------------

void Layout::block_indices 
(int ib, int * ibx, int * iby, int * ibz) const throw()
{
  if (ibx) (*ibx) = ib % block_count_[0];
  if (iby) (*iby) = (ib / block_count_[0]) % block_count_[1];
  if (ibz) (*ibz) = ib / (block_count_[0]*block_count_[1]);
}

//----------------------------------------------------------------------

void Layout::block_indices (int ip, int ibl, int * ibx, int * iby, int * ibz) const throw()
{
  int ibg = global_index(ip,ibl);
  block_indices(ibg,ibx,iby,ibz);
}
