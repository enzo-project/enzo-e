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
    p | nix_;
    p | niy_;
    p | niz_;
    p | field_group_;
    p | num_fields_;
    p | field_values_;
    p | volume_ratio_;
    p | spheres_;
  }

  virtual ~EnzoLevelArray()
  {
    CkPrintf ("~EnzoLevelArray() %d %d %d\n",
              thisIndex[0],thisIndex[1],thisIndex[2]);
    fflush(stdout);
  }
  
  /// Send data request to containing EnzoBlock in level_base
  void p_request_data ();

  /// Accept requested data from blocks
  void p_transfer_data (Index, int nf, enzo_float * field_data_list );

  /// Apply the inference method on the arrays, synchronizing
  /// with EnzoSimulation[0] afterwards
  void apply_inference();

protected: // functions

  /// Return the index limits of octree root blocks intersecting this
  /// level array
  void intersecting_root_blocks_
  (Block * block, int im3[3], int ip3[3], int array_size[3]) const;

  void coarsen_
  (enzo_float * ac, int mcx, int mcy, int mcz, int ncx, int ncy, int ncz,
   const enzo_float * af, int mfx, int mfy, int mfz, int nfx, int nfy, int nfz);

  void interpolate_
  (enzo_float * af,
   int mfx, int mfy, int mfz, int nfx, int nfy, int nfz, int efx, int efy, int efz,
   const enzo_float * ac,
   int mcx, int mcy, int mcz, int ncx, int ncy, int ncz, int ecx, int ecy, int ecz);

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

  /// Dimensions of the inference arrays
  int nix_, niy_, niz_;

  /// Field group defining fields stored in level arrays
  std::string field_group_;

  /// Number of fields in the field_group_
  int num_fields_;

  /// Arrays of fields
  std::vector< std::vector < enzo_float > > field_values_;

  /// Variable for keeping track of volume of incoming data
  float volume_ratio_;

  /// List of spheres
  std::vector<ObjectSphere> spheres_;
};

#endif /* ENZO_IO_ENZO_LEVEL_ARRAY_HPP */

