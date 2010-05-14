// $Id: parallel_layout.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef PARALLEL_LAYOUT_HPP
#define PARALLEL_LAYOUT_HPP

/// @file     parallel_layout.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the Layout class

#include <vector>

class Layout {

  /// @class    Layout
  /// @ingroup  Layout
  /// @brief Specify how a Patch is partitioned into process, thread, and computational blacks

public: // interface

  /// (*) Create a Layout of the given dimensionality, defaulting to serial
  Layout() throw();

  /// ( ) Destructor
  //  ~Layout() throw();

  /// ( ) Copy constructor
  //  Layout(const Layout & layout) throw() throw();

  /// ( ) Assignment operator
  //  Layout & operator= (const Layout & layout) throw() throw();

  /// Set how many process blocks the Layout will be partitioned into 
  void set_process_blocks(int p0, int p1, int p2) throw()
  { process_blocks_[0] = p0;
    process_blocks_[1] = p1;
    process_blocks_[2] = p2; };
    

  /// Set how thread blocks each process block is partitioned into 
  void set_thread_blocks(int t0, int t1, int t2) throw()
  { thread_blocks_[0] = t0;
    thread_blocks_[1] = t1;
    thread_blocks_[2] = t2; };

  /// Set how many compute blocks each thread block is partitioned into  
  void set_data_blocks(int d0, int d1, int d2) throw()
  { data_blocks_[0] = d0;
    data_blocks_[1] = d1;
    data_blocks_[2] = d2; };

  /// Return the number of processes in the Layout 
  int process_count () throw()
  { return 
      process_blocks_[0]*
      process_blocks_[1]*
      process_blocks_[2]; };

  /// Return the number of threads per process in the Layout 
  int thread_count () throw()
  { return 
      thread_blocks_[0]*
      thread_blocks_[1]*
      thread_blocks_[2]; };

  /// Return the number of data blocks per thread in the Layout  
  int data_blocks_per_thread () throw()
  {
    return 
      data_blocks_[0]*
      data_blocks_[1]*
      data_blocks_[2];
  }

  /// Return the number of data blocks per process in the Layout  
  int data_blocks_per_process () throw()
  {
    return thread_count() * data_blocks_per_thread();
  };

  /// Return the relative process block in the given direction from
  //  the given (process,thread,data) block
  int neighbor_process_block (int process, int thread, int block,
			      int axis, int dir)  throw();

  /// Return the relative thead block in the given direction from the
  /// given (process,thread,data) block
  int neighbor_thread_block (int process, int thread, int block,
			     int axis, int dir) throw(); 

  /// Return the relative data block in the given direction from the
  /// given (process,thread,data) block
  int neighbor_data_block (int process, int thread, int block,
			   int axis, int dir) throw(); 

  /// Return the bounds associated with the given
  /// (process,thread,block) relative to the Patch block 0 < x,y,z < 1
  void box_extent (int process, int thread, int block,
		   double lower_extent[3],    
		   double upper_extent[3]) ;

  /// Return the index range assigned to the given
  /// (process,thread,block) given a Patch array size
  void array_indices (int process, int thread, int block,
		      int lower_index[3],
		      int upper_index[3]) throw ();

private: // attributes

  ///  number of distributed memory processes
  int process_blocks_[3];

  /// number of shared memory threads per process
  int thread_blocks_[3];

  /// number of compute blocks per thread
  int data_blocks_[3];

};
#endif /* PARALLEL_LAYOUT_HPP */
