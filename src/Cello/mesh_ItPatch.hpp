// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItPatch.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Mesh] Declaration of the ItPatch iterator

#ifndef MESH_IT_PATCH_HPP
#define MESH_IT_PATCH_HPP

class ItPatch : public It<Patch> {

  /// @class    ItPatch
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator over Patches in a Hierarchy

public: // interface

  /// Create an ItPatch object
  ItPatch (const Hierarchy * hierarchy) throw ();

  /// Delete the ItPatch object
  virtual ~ItPatch () throw ();
  
  /// Iterate through all local Patches in the Hierarchy
  Patch * operator++ () throw();

  /// Return whether the iteration is complete
  bool done() const throw();

private: // attributes

  /// The Hierarchy being iterated over
  const Hierarchy * hierarchy_;

  /// Index of the current local Patch plus 1, or 0 if between iterations
  /// Always in the range 0 <= index1_ <= number of local patches
  size_t index1_;
};

#endif /* MESH_IT_PATCH_HPP */
