// $Id: array_layout.hpp 1258 2010-03-02 01:07:36Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ARRAY_LAYOUT_HPP
#define ARRAY_LAYOUT_HPP

/// @file     array_layout.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Declaration of the Layout class

#include <vector>

class Layout {

  /// @class    Layout
  /// @ingroup  Array
  /// @brief Specify how the Array Block s are distributed across
  /// processes and threads

public: // interface

  /// (*) Create a Layout of the given dimensionality, defaulting to serial
  Layout(int dimension);

  /// ( ) Destructor
  ~Layout();

  /// ( ) Copy constructor
  Layout(const Layout & layout) throw();

  /// ( ) Assignment operator
  Layout & operator= (const Layout & layout) throw();

  //----------------------------------------------------------------------
  // Array level 
  //----------------------------------------------------------------------

  /// (*) Set the size of arrays (dimension used only for checking
  /// array bounds)
  void set_array ( int dimension, int array_size[] );

  /// (*) Return the size of arrays (dimension used only for checking
  /// array bounds)
  void array_size ( int dimension, int array_size[] );

  //----------------------------------------------------------------------
  // Physical processes
  //----------------------------------------------------------------------

  /// (*) Set first physical process and number of physical processes
  void set_processes(int process_first, int process_count)
  { process_first_ = process_first;  process_count_ = process_count;  };

  /// (*) Return the first physical process assigned to this layout
  int process_first () { return process_first_; };

  /// (*) Return the number of physical processes assigned to this layout
  int process_count () { return process_count_; };

  //----------------------------------------------------------------------
  // Physical threads
  //----------------------------------------------------------------------

  /// (*) Set first physical thread and number of physical threads
  void set_threads(int thread_first, int thread_count)
  { thread_first_ = thread_first; thread_count_ = thread_count;  };

  /// (*) Return the first physical thread assigned to this layout
  int thread_first() { return thread_first_; };

  /// (*) Return the number of physical threads assigned to this layout
  int thread_count() { return thread_count_; };

  //----------------------------------------------------------------------
  // Process block level
  //----------------------------------------------------------------------

  /// (*) Set virtual process partitioning of the array
  void set_process_blocks (int dimension, int process_block_count[]);

  /// (*) Return the number of local process blocks
  int process_block_count (int process_rank) const;

  /// ( ) Return the index of the local process block for the given block index
  void process_block_indices
  (
   int dimension, 
   int process_rank, 
   int process_block_number, 
   int process_block_indices[] );

  /// ( ) Return the process rank and process block number for the given
  /// process block indices (invers of process_block_indices())
  void process_block_number
  (
   int dimension, 
   int * process_rank, 
   int * process_block_number, 
   int process_block_indices[] );

  //----------------------------------------------------------------------
  // Thread block level
  //----------------------------------------------------------------------

  /// ( ) Set virtual thread partitioning of a process block
  void set_thread_blocks (int dimension, int thread_block_count[]);

  /// ( ) Return the number of thread blocks
  int thread_block_count ( int index_process_block ) const;

  /// ( ) Return the index of the thread block for the given thread block index
  int thread_block_indices
  (
   int index_process_block, 
   int index_thread_block, 
   int dimension, 
   int thread_block_indices[]
   );

private: // attributes

  ///  dimensionality of the `Array`  
  int dimension_;

  ///  first physical process 
  int process_first_;

  ///  last physical process 
  int process_count_;

  ///  first physical thread 
  int thread_first_;

  ///  last physical thread 
  int thread_count_;

  ///  Size of the entire array
  int * array_size_;

  /// Product of process_blocks_ elements
  int process_block_count_;

  /// Process partitioning of array
  int * process_blocks_;

  /// Product of thread_blocks_ elements
  int thread_block_count_;

  /// Thread partitioning of process blocks
  int * thread_blocks_;

};
#endif /* ARRAY_LAYOUT_HPP */
