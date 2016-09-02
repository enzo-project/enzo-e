// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MappingArray.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-04-22
/// @brief    [\ref Parallel] Declaration of the MappingArray class

#ifndef CHARM_MAPPING_ARRAY_HPP
#define CHARM_MAPPING_ARRAY_HPP

#include "simulation.decl.h"

class MappingArray: public CkArrayMap {

  /// @class    MappingArray
  /// @ingroup  Charm
  /// @brief    [\ref Parallel] Class for mapping Blocks to processors
  ///
  /// This class defines how to map a 3D array of Charm++ chares to
  /// processes

public:

  MappingArray(int nx, int ny, int nz);

  int procNum(int, const CkArrayIndex &idx);

  /// CHARM++ migration constructor for PUP::able
  MappingArray (CkMigrateMessage *m)
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

#endif /* CHARM_MAPPING_ARRAY_HPP */

