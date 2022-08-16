// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoLevelArray.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2022-07-29
/// @brief    [\ref Io] Declaration of the EnzoLevelArray class

#ifndef ENZO_IO_ENZO_LEVEL_ARRAY_HPP
#define ENZO_IO_ENZO_LEVEL_ARRAY_HPP

class EnzoLevelArray : public CBase_EnzoLevelArray {

  /// @class    EnzoLevelArray
  /// @ingroup  Enzo
  /// @brief    [\ref Io] 

public: // interface

  /// Constructors
  EnzoLevelArray()
  {
    CkPrintf ("TRACE EnzoLevelArray %d %d %d\n",thisIndex.x,thisIndex.y,thisIndex.z);
  }
  
  /// CHARM++ migration constructor
  EnzoLevelArray(CkMigrateMessage *m) : CBase_EnzoLevelArray(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP;
  }


protected: // functions

  /// Return the index limits of octree root blocks intersecting this
  /// level array
  void intersecting_root_blocks_(Block * block, int im3[3], int ip3[3], int array_size[3]) const;

private: // attributes

};

#endif /* ENZO_IO_ENZO_LEVEL_ARRAY_HPP */

