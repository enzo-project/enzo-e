// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_OrderingHilbert.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-02-01
/// @brief    [\ref Mesh] Declaration of the OrderingHilbert class

#ifndef MESH_ORDERING_HILBERT_HPP
#define MESH_ORDERING_HILBERT_HPP

class OrderingHilbert : public Ordering {

  /// @class    OrderingHilbert
  /// @ingroup  Mesh
  /// @brief    [\ref Mesh] 

public: // interface

  /// Constructor
  OrderingHilbert() throw()
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(OrderingHilbert);

  /// CHARM++ migration constructor
  OrderingHilbert(CkMigrateMessage *m) : Ordering(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP; Ordering::pup(p); }

public:  // virtual methods

  /// Return the index of the given Block
  virtual int index (Block * block) const
  { return -1; }

  /// Return the number of indices in the OrderingRootBlocks
  virtual int num_indices () const
  { return 0; }

private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* MESH_ORDERING_HILBERT_HPP */

