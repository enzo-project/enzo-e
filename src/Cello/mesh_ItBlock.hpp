// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_ItBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Feb  1 16:46:01 PST 2011
/// @brief    [\ref Mesh] Declaration of the ItBlock iterator

#ifndef MESH_IT_BLOCK_HPP
#define MESH_IT_BLOCK_HPP

class ItBlock : public It<CommBlock> {

  /// @class    ItBlock
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] Iterator over local CommBlocks in a Patch

public: // interface

  /// Create an ItBlock object
#ifdef REMOVE_PATCH
  ItBlock (const Hierarchy * hierarchy) throw ();
#else /* REMOVE_PATCH */
  ItBlock (const Patch * patch) throw ();
#endif /* REMOVE_PATCH */


  /// Delete the ItBlock object
  virtual ~ItBlock () throw ();

#ifdef CONFIG_USE_CHARM
  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
    It<CommBlock>::pup(p);
    WARNING("ItBlock::pup","skipping patch_ (aliased pointer)");
    // p | *patch_;
  }
#endif
  
  /// Iterate through all local CommBlocks in the Patch
  CommBlock * operator++ () throw();

  /// Return whether the iteration is complete
  bool done() const throw();

private: // attributes

#ifdef REMOVE_PATCH
  Hierarchy * hierarchy_;
#else /* REMOVE_PATCH */
  Patch * patch_;
#endif /* REMOVE_PATCH */


};

#endif /* MESH_IT_BLOCK_HPP */
