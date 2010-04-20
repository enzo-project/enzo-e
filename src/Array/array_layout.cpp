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
  
  array_size_     = new int [dimension_];
  process_blocks_ = new int [dimension_];
  thread_blocks_  = new int [dimension_];

  for (int i=0; i<dimension; i++) {
    array_size_[i]     = 1;
    process_blocks_[i] = 1;
    thread_blocks_[i]  = 1;
  }
  Parallel * parallel = Parallel::instance();
  process_rank_ = parallel->process_rank();
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
 int process_blocks[]
 )
{
  int i;
  process_block_count_ = 1;
  for (i=0; i<MIN(dimension,dimension_); i++) {
    process_blocks_[i] = process_blocks[i];
    process_block_count_ *= process_blocks[i];
  }
  // Fill any extra elements of output process_blocks[] with 1's
  for (i=MIN(dimension,dimension_); i<dimension_; i++) {
    process_blocks_[i] = 1;
  }
 
}

//----------------------------------------------------------------------

int Layout::process_block_count () const
/// @todo store local process block count result
/// @todo store global process block count (nvp)
{

  int process_rank_local = process_rank_ - process_first_;

  int block_count;

  if (0 <= process_rank_local && 
      process_rank_local < process_count_) {

    // compute process block range on this process

    int block_index_first 
      = process_rank_local     * process_block_count_ / process_count_;

    int block_index_last  
      = (process_rank_local + 1) * process_block_count_ / process_count_;

    block_count = block_index_last - block_index_first;
    
  } else {
    // Array does not belong on this process
    block_count = 0;
  }

  return block_count;
}

//----------------------------------------------------------------------

void Layout::process_block_indices
(
 int block_offset, 
 int dimension, 
 int process_block_index[] )
{
  ASSERT ("Layout::process_block_indices",
	  "block_offset out of range",
	  0 <= block_offset &&
	  block_offset < process_block_count_);

  int process_rank_local = process_rank_ - process_first_;

  if (0 <= process_rank_local && 
      process_rank_local < process_count_) {

    // compute process block range on this process

    int block_index_first 
      = process_rank_local * process_block_count_ / process_count_;

    int block_index = block_index_first + block_offset;

    for (int i=0; i<MIN(dimension,dimension_); i++) {
      int index = block_index % process_blocks_[i];
      block_index = (block_index - index) / process_blocks_[i];
      process_block_index[i] = index;
    }
    // pad any extra indices with 0's
    for (int i=MIN(dimension,dimension_)+1; i<dimension; i++) {
      process_block_index[i] = 0;
    }
  } else {
    WARNING_MESSAGE("process_block_indices","Process out of range--ignoring");
  }
}

//----------------------------------------------------------------------

void Layout::set_thread_blocks (int dimension, int thread_blocks[])
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
