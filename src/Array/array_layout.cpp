// $Id: array_layout.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     array_layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the Layout class

//----------------------------------------------------------------------

#include "error.hpp"
#include "parallel.hpp"
#include "array_layout.hpp"

//----------------------------------------------------------------------

Layout::Layout(int dimension)
  : dimension_(dimension),
    process_first_(0),
    process_count_(1),
    thread_first_(0),
    thread_count_(1)
    
{

  ASSERT ("Layout::Layout",
	  "No funny dimensions allowed",
	  1 <= dimension && dimension <= 3);
  
  array_size_          = new int [dimension_];
  process_block_count_ = new int [dimension_];
  thread_block_count_  = new int [dimension_];

  for (int i=0; i<dimension; i++) {
    array_size_[i]         = 1;
    process_block_count_[i] = 1;
    thread_block_count_[i]  = 1;
  }
}

//----------------------------------------------------------------------

Layout::~Layout()
{
}

//----------------------------------------------------------------------

void Layout::set_array
(
 int dimension, 
 int array_size[]
 )
{
  if (dimension > dimension_) {
    WARNING_MESSAGE("Layout::set_array","dimension > dimension_");
  }
  int i;
  for (i=0; i<MIN(dimension,dimension_); i++) {
    array_size_[i] = array_size[i];
  }
  // Fill any extra elements of Layout's array size with 1's
  for (i=MIN(dimension,dimension_); i<dimension_; i++) {
    array_size_[i] = 1;
  }
}

//----------------------------------------------------------------------

void Layout::array_size 
(
 int dimension, 
 int array_size[]
 )
{
  if (dimension > dimension_) {
    WARNING_MESSAGE("Layout::array_size","dimension > dimension_");
  }
  int i;
  for (i=0; i<MIN(dimension,dimension_); i++) {
    array_size[i] = array_size_[i];
  }
  // Fill any extra elements of output array_size[] with 1's
  for (i=MIN(dimension,dimension_); i<dimension; i++) {
    array_size[i] = 1;
  }
}

//----------------------------------------------------------------------

/// Set the range process_first to process_max+1 of processors
void Layout::set_processes
(
 int process_first, 
 int process_count
 )
{
  process_first_ = process_first;
  process_count_ = process_count;
}

//----------------------------------------------------------------------

/// Set the range thread_first to thread_max+1 of threads
void Layout::set_threads
(
 int thread_first, 
 int thread_count
 )
{
  thread_first_ = thread_first;
  thread_count_ = thread_count;
}

//----------------------------------------------------------------------

void Layout::set_process_blocks 
(
 int dimension, 
 int process_block_count[]
 )
{
  int i;
  for (i=0; i<MIN(dimension,dimension_); i++) {
    process_block_count_[i] = process_block_count[i];
  }
  // Fill any extra elements of output process_block_count[] with 1's
  for (i=MIN(dimension,dimension_); i<dimension_; i++) {
    process_block_count_[i] = 1;
  }
 
}

//----------------------------------------------------------------------

int Layout::process_block_count () const
{
  // threads per process
  Parallel * parallel = Parallel::instance();
  int ip = parallel->process_rank();

  int nb = 0;
  if (process_first_ <= ip && ip < process_first_ + process_count_) {
    // array size             array_size_[]
    // virtual process counts process_block_count_[]
    // physical process       process_count_
    //
    // vp = proces
    // In range
    // int ip0 = ip - process_first_;
    // float ia0 = ip0 * process_count_ / array_size_[0];
    int nvp = 1;
    for (int i=0; i<dimension_; i++) {
      nvp *= process_block_count_[i];
    }
    int ip_first = ip*nvp / process_count_;
    int ip_last  = (ip+1)*nvp / process_count_;
    printf ("vp_first = %d  vp_last = %d\n",ip_first,ip_last);
    nb = ip_last - ip_first;
    
  } else {
    // Out of range
    nb = 0;
  }

  printf ("ip nb = %d %d\n",ip,nb);
  return nb;
}

//----------------------------------------------------------------------

int Layout::process_block_indices
(
 int index_process_block, 
 int dimension, 
 int process_block_indices[] )
{
  INCOMPLETE_MESSAGE("Layout::process_block_indices","");
  return 0;
}

//----------------------------------------------------------------------

void Layout::set_thread_blocks (int dimension, int thread_block_count[])
{
  INCOMPLETE_MESSAGE("Layout::set_thread_blocks","");
}

//----------------------------------------------------------------------

int Layout::thread_block_count ( int index_process_block ) const
{
  INCOMPLETE_MESSAGE("Layout::thread_block_count","");
  return 0;
}

//----------------------------------------------------------------------

int Layout::thread_block_indices
(
 int index_process_block, 
 int index_thread_block, 
 int dimension, 
 int thread_block_indices[]
 )
{
  INCOMPLETE_MESSAGE("Layout::thread_block_indices","");
  return 0;
}

//----------------------------------------------------------------------
