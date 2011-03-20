// $Id: mesh_ItBlock.hpp 1942 2011-01-20 00:53:45Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_IT_BLOCKS_HPP
#define MESH_IT_BLOCKS_HPP

/// @file     mesh_ItBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Move creation of iterator to iterated object: Mesh::create_iter() (factor method)
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Mesh] Declaration of the ItBlock iterator

class ItBlock {

  /// @class    ItBlock
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator over Blocks in a Patch

public: // interface

  /// Create an ItBlock object
  ItBlock (Patch * patch) throw ();

  /// Delete the ItBlock object
  ~ItBlock () throw ();
  
  /// Iterate through all local Blocks in the Patch
  Block * operator++ () throw();

  /// Return whether the iteration is complete
  bool done() const throw();

private: // attributes

  /// The Patch being iterated over
  Patch * patch_;

  /// Index of the current local Block plus 1, or 0 if between iterations
  /// Always in the range 0 <= index1_ <= number of local blocks
  size_t index1_;
};

#endif /* MESH_IT_BLOCKS_HPP */
