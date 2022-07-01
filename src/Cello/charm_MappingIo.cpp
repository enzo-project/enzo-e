// See LICENSE_CELLO file for license and copyright information

/// @file     charm_MappingIo.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-04-12
/// @brief    Mapping of Charm++ array Index to processors

#include "charm.hpp"

//======================================================================

MappingIo::MappingIo(int count)
  :  CkArrayMap()
{
  count_ = count;
}

//----------------------------------------------------------------------

int MappingIo::procNum(int, const CkArrayIndex &idx) {

  int index_block = *(int*)idx.data();
#ifdef CONFIG_SMP_MODE
  return count_? (long long) index_block*CkNumNodes()/count_ : index_block;
#else
  return count_? (long long) index_block*CkNumPes()/count_ : index_block;
#endif
}

