// $Id: mesh_ItBlocks.hpp 1942 2011-01-20 00:53:45Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_IT_BLOCKS_HPP
#define MESH_IT_BLOCKS_HPP

/// @file     mesh_ItBlocks.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Implement operator++() to iterate over all Patches not just root
/// @todo     Implement iterators using templates to avoid void *
/// @todo     Move ItBlocks to Mesh component
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Mesh] Declaration of the ItBlocks iterator

class ItBlocks {

  /// @class    ItBlocks
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator over Blocks in a Patch

public: // interface

  /// Create an ItBlocks object
  ItBlocks (Patch * patch) throw ();

  /// Delete the ItBlocks object
  ~ItBlocks () throw ();
  
  /// Iterate through all Patches in the patch (currently only root!)
  DataBlock * operator++ ();

private: // attributes

  /// The Patch being iterated over
  Patch * patch_;

  /// Index to the current DataBlock
  size_t curr_;
};

#endif /* MESH_IT_BLOCKS_HPP */
