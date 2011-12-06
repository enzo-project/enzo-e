// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_BlockReduce.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-08-10
/// @todo     Change BlockReduce from Block->Mesh to Block->Patch reductions
/// @todo     Add PatchReduce for Patch->Mesh reductions
/// @brief    [\ref Mesh] Declaration of the BlockReduce class

#ifndef MESH_BLOCK_REDUCE_HPP
#define MESH_BLOCK_REDUCE_HPP

#ifdef CONFIG_USE_CHARM

class BlockReduce {

  /// @class    BlockReduce
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] CHARM chare class for parallel reductions from Blocks to Mesh level

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

  /// Prepare counters
  int count_prepare_;

  /// Output counters
  int count_output_;

  /// Prepare reduction variables
  double dt_hierarchy_;
  int stop_hierarchy_;

};

#endif /* CONFIG_USE_CHARM */

#endif /* MESH_BLOCK_REDUCE_HPP */

