// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MappingTree.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-22
/// @brief    [\ref Parallel] Declaration of the MappingTree class

#ifndef CHARM_MAPPING_TREE_HPP
#define CHARM_MAPPING_TREE_HPP

#include "simulation.decl.h"

class MappingTree: public CkArrayMap {

  /// @class    MappingTree
  /// @ingroup  Charm
  /// @brief    [\ref Parallel] Class for mapping Blocks to processors
  ///
  /// This class defines how to map a 3D array of Charm++ chares to
  /// processes

public:

  MappingTree(int nx, int ny, int nz);

  int procNum(int, const CkArrayIndex &idx);

  /// CHARM++ migration constructor for PUP::able
  MappingTree (CkMigrateMessage *m)
    : CkArrayMap(m),
      nx_(0),ny_(0),nz_(0)
  { }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    CkArrayMap::pup(p);
    // NOTE: change this function whenever attributes change
    p | nx_;
    p | ny_;
    p | nz_;
  }

private:

  int nx_, ny_, nz_;

};

#endif /* CHARM_MAPPING_TREE_HPP */

