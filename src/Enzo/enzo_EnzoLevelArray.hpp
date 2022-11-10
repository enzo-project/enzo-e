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
  EnzoLevelArray(std::string field_group,
                 int level_base, int level_array, int level_infer,
                 int nax, int nay=1, int naz=1);
  
  /// CHARM++ migration constructor
  EnzoLevelArray(CkMigrateMessage *m) : CBase_EnzoLevelArray(m) {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p) 
  { TRACEPUP;
    p | level_base_;
    p | level_array_;
    p | level_infer_;
    p | nax_;
    p | nay_;
    p | naz_;
    p | field_group_;
    p | field_values_;
    p | volume_ratio_;
  }

  /// Send data request to containing EnzoBlock in level_base
  void p_request_data ();

  /// Accept requested data from blocks
  void p_transfer_data (Index, int nf, enzo_float * field_data_list );


protected: // functions

  /// Return the index limits of octree root blocks intersecting this
  /// level array
  void intersecting_root_blocks_
  (Block * block, int im3[3], int ip3[3], int array_size[3]) const;

  void coarsen_
  (enzo_float * ac, int mcx, int mcy, int mcz,
   enzo_float * af, int mfx, int mfy, int mfz);

  void interpolate_
  (enzo_float * ac, int mcx, int mcy, int mcz,
   enzo_float * af, int mfx, int mfy, int mfz);

private: // attributes

  /// AMR level of blocks associated with this array
  /// 0 <= level_base_ <= level_array_
  int level_array_;

  /// AMR level of blocks used as base for communicating with blocks
  /// 0 <= level_base_ <= level_array_
  int level_base_;

  /// AMR level where AMR blocks have the same resolution as the
  /// inference arrays
  int level_infer_;

  /// Dimensions of this array
  int nax_, nay_, naz_;

  /// Field group defining fields stored in level arrays
  std::string field_group_;

  /// Arrays of fields
  std::vector< std::vector < enzo_float > > field_values_;

  /// Variable for keeping track of volomue of incoming data
  float volume_ratio_;
};

#endif /* ENZO_IO_ENZO_LEVEL_ARRAY_HPP */

