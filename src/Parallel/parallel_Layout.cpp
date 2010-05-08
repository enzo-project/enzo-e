// $Id: parallel_layout.cpp 1262 2010-03-03 15:44:05Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     parallel_layout.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the Layout class

//----------------------------------------------------------------------

#include "error.hpp"
#include "parallel.hpp"

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
}

//----------------------------------------------------------------------

Layout::~Layout()
{
  delete [] array_size_;
  delete [] process_blocks_;
  delete [] thread_blocks_;
}

//----------------------------------------------------------------------

Layout::Layout(const Layout & layout) throw()
  : dimension_(layout.dimension_),
    process_first_(layout.process_first_),
    process_count_(layout.process_count_),
    thread_first_(layout.thread_first_),
    thread_count_(layout.thread_count_)
{
  array_size_     = new int [dimension_];
  process_blocks_ = new int [dimension_];
  thread_blocks_  = new int [dimension_];

  for (int i=0; i<dimension_; i++) {
    array_size_[i]     = layout.array_size_[i];
    process_blocks_[i] = layout.process_blocks_[i];
    thread_blocks_[i]  = layout.thread_blocks_[i];
  }
}

//----------------------------------------------------------------------

Layout & Layout::operator= (const Layout & layout) throw()
{
  dimension_     = layout.dimension_;
  process_first_ = layout.process_first_;
  process_count_ = layout.process_count_;
  thread_first_  = layout.thread_first_;
  thread_count_  = layout.thread_count_;

  array_size_     = new int [dimension_];
  process_blocks_ = new int [dimension_];
  thread_blocks_  = new int [dimension_];

  for (int i=0; i<dimension_; i++) {
    array_size_[i]     = layout.array_size_[i];
    process_blocks_[i] = layout.process_blocks_[i];
    thread_blocks_[i]  = layout.thread_blocks_[i];
  }

  return *this;
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

void Layout::set_process_blocks 
(
 int dimension, 
 int process_blocks[]
 )
{
  int i;
  for (i=0; i<MIN(dimension,dimension_); i++) {
    process_blocks_[i] = process_blocks[i];
  }
  // Fill any extra elements of output process_blocks[] with 1's
  for (i=MIN(dimension,dimension_); i<dimension_; i++) {
    process_blocks_[i] = 1;
  }
 
}

//----------------------------------------------------------------------

int Layout::process_block_count (int process_rank) const
/// @todo store local process block count result
/// @todo store global process block count (nvp)
{

  int process_rank_local = process_rank - process_first_;

  int block_count = 0;

  if (0 <= process_rank_local && 
      process_rank_local < process_count_) {

    // compute process block range on this process

    int process_block_count = 1;
    for (int i=0; i<dimension_; i++) process_block_count *= process_blocks_[i];

    int block_index_first 
      = process_rank_local     * process_block_count / process_count_;

    int block_index_last  
      = (process_rank_local + 1) * process_block_count / process_count_;

    block_count = block_index_last - block_index_first;
    
  }

  return block_count;
}

//----------------------------------------------------------------------

void Layout::process_block_indices
(
 int dimension, 
 int process_rank,
 int process_block_number, 
 int process_block_index[] )
{
  int process_block_count = 1;
  for (int i=0; i<dimension; i++) process_block_count *= process_blocks_[i];
  ASSERT ("Layout::process_block_indices",
	  "process_block_number out of range",
	  0 <= process_block_number &&
	  process_block_number < process_block_count);

  int process_rank_local = process_rank - process_first_;

  if (0 <= process_rank_local && 
      process_rank_local < process_count_) {

    // compute process block range on this process

    int block_index_first 
      = process_rank_local * process_block_count / process_count_;

    int block_index = block_index_first + process_block_number;

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

void Layout::process_block_number
(
 int dimension, 
 int * process_rank,
 int * process_block_number, 
 int process_block_index[] )
{
  int d = MIN(dimension,dimension_)-1;
  int block_index = process_block_index[d];
  for (int i=d-1; i>=0; i--) {
    block_index = block_index * process_blocks_[i] + process_block_index[i];
  }

  // compute process_block_count
  int process_block_count = 1;
  for (int i=0; i<dimension; i++) process_block_count *= process_blocks_[i];

  // compute *process_rank

  int process_rank_local
    = block_index * process_count_ / process_block_count;

  int block_index_first =
    (process_rank_local + 1) * process_block_count / process_count_;

  while ( block_index >= block_index_first) {
    ++ process_rank_local;
    block_index_first =
      (process_rank_local + 1) * process_block_count / process_count_;
  }

  *process_rank         = process_rank_local + process_first_;

  block_index_first =
    process_rank_local * process_block_count / process_count_;
  *process_block_number = block_index - block_index_first;

  printf ("%d %d %d  %d   %d  %d\n",
	  process_block_index[0],
	  process_block_index[1],
	  process_block_index[2],
	  block_index,
	  *process_rank,
	  *process_block_number);
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
