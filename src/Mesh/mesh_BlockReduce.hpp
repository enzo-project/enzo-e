// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_BLOCK_REDUCE_HPP
#define MESH_BLOCK_REDUCE_HPP

/// @file     mesh_BlockReduce.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    [\ref Mesh] Brief description of file mesh_BlockReduce.hpp

#ifdef CONFIG_USE_CHARM

class BlockReduce {

  /// @class    BlockReduce
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Brief description of class BlockReduce.

public: // interface

  /// Constructor
  BlockReduce();

  /// Reduce block dt and stop, then proceed with cycle
  void p_prepare(int count, int cycle, double time,
		 double dt_block, int stop_block);

  /// Reduce output from simulation
  void p_output_reduce(int count);

private: // functions


private: // attributes

  /// Output counters
  int count_output_;

  /// Prepare counters
  int count_prepare_;

  /// Prepare reduction variables
  double dt_hierarchy_;
  int stop_hierarchy_;

};

#endif /* CONFIG_USE_CHARM */

#endif /* MESH_BLOCK_REDUCE_HPP */

