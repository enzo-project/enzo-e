// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Mesh] Declaration of the ItBlock iterator

#ifndef MESH_IT_BLOCK_HPP
#define MESH_IT_BLOCK_HPP

class ItBlock : public It<Block> {

  /// @class    ItBlock
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator over local Blocks in a Patch

public: // interface

  /// Create an ItBlock object
  ItBlock (const Patch * patch) throw ();

  /// Delete the ItBlock object
  virtual ~ItBlock () throw ();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    It<Block>::pup(p);
    p | *patch_;
  }
#endif
  
  /// Iterate through all local Blocks in the Patch
  Block * operator++ () throw();

  /// Return whether the iteration is complete
  bool done() const throw();

private: // attributes

  /// The Patch being iterated over
  Patch * patch_;

};

#endif /* MESH_IT_BLOCK_HPP */
