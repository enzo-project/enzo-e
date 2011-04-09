// $Id: enzo_EnzoBlock.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

#ifndef ENZO_ENZO_BLOCK_HPP
#define ENZO_ENZO_BLOCK_HPP

/// @file     enzo_EnzoBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:01:51 PST 2011
/// @todo     Change public attributes to private
/// @todo     Dynamically allocate arrays
/// @brief    [\ref Enzo] Declaration of the EnzoBlock class

//----------------------------------------------------------------------

class EnzoBlock : public Block , public CBase_EnzoBlock  
{

  /// @class    EnzoBlock
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] An EnzoBlock is a Block with Enzo data

public: // interface

  /// Initialize the BlockCharm chare
  EnzoBlock
  ( int nx, int ny, int nz,
    double xm, double ym, double zm,
    double hx, double hy, double hz,
    int num_field_blocks) throw();

  /// Initialize a migrated Block
  EnzoBlock (CkMigrateMessage *m) {CkPrintf("Oops\n");};

  /// Initialize an empty Block
  EnzoBlock() {CkPrintf("Oops\n");};
};

#endif /* ENZO_ENZO_BLOCK_HPP */

