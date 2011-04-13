// $Id: mesh_ItBlockLocal.hpp 1942 2011-01-20 00:53:45Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef MESH_IT_BLOCKS_HPP
#define MESH_IT_BLOCKS_HPP

/// @file     mesh_ItBlockLocal.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @todo     Move creation of iterator to iterated object: Mesh::create_iter() (factor method)
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Mesh] Declaration of the ItBlockLocal iterator

#ifndef CONFIG_USE_CHARM

class ItBlockLocal {

  /// @class    ItBlockLocal
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator over local Blocks! in a Patch

public: // interface

  /// Create an ItBlockLocal object
  ItBlockLocal (Patch * patch) throw ();

  /// Delete the ItBlockLocal object
  ~ItBlockLocal () throw ();
  
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

#endif

#endif /* MESH_IT_BLOCKS_HPP */
