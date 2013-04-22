// See LICENSE_CELLO file for license and copyright information

/// @file     charm_ArrayMap.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Charm] Declaration of the ArrayMap class
///

#ifndef CHARM_ARRAY_MAP_HPP
#define CHARM_ARRAY_MAP_HPP

#include "simulation.decl.h"

class ArrayMap: public CBase_ArrayMap {

  /// @class    ArrayMap
  /// @ingroup  Charm
  /// @brief    [\ref Charm] 

public:
  int *mapping;

  ArrayMap(int nx, int ny, int nz);
  ~ArrayMap();
  int procNum(int, const CkArrayIndex &idx);

private:

  int nx_, ny_, nz_;

};

#endif /* CHARM_ARRAY_MAP_HPP */

