// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MappingIo.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-04-12
/// @brief    [\ref Parallel] Declaration of the MappingIo class

#ifndef CHARM_MAPPING_IO_HPP
#define CHARM_MAPPING_IO_HPP

#include "simulation.decl.h"

class MappingIo: public CkArrayMap {

  /// @class    MappingIo
  /// @ingroup  Charm
  /// @brief    [\ref Parallel] Class for mapping IoWriter and IoReader elements to processors

public:

  MappingIo(int count);

  int procNum(int, const CkArrayIndex &idx);

  /// CHARM++ migration constructor for PUP::able
  MappingIo (CkMigrateMessage *m)
    : CkArrayMap(m),
      count_(0)
  { }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    CkArrayMap::pup(p);
    // NOTE: change this function whenever attributes change
    p | count_;
  }

private:

  int count_;

};

#endif /* CHARM_MAPPING_IO_HPP */

