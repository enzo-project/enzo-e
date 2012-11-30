// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCommBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-28
/// @brief    [\ref Enzo] Declaration of the EnzoCommBlock class

#ifndef ENZO_ENZO_COMM_BLOCK_HPP
#define ENZO_ENZO_COMM_BLOCK_HPP

class EnzoCommBlock : CommBlock {

  /// @class    EnzoCommBlock
  /// @ingroup  Comm
  /// @brief    [\ref Comm] 

public: // interface

  /// Constructor
  EnzoCommBlock
  (
   int ibx, int iby, int ibz,
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xmp, double ymp, double zmp,
   double xb, double yb, double zb,
   CProxy_Patch proxy_patch,
   int patch_id,
   int patch_rank,
   int num_field_blocks
   ) throw();

  EnzoCommBlock
  (
   int nbx, int nby, int nbz,
   int nx, int ny, int nz,
   double xmp, double ymp, double zmp,
   double xb, double yb, double zb,
   CProxy_Patch proxy_patch,
   int patch_id,
   int patch_rank,
   int num_field_blocks
) throw();

  /// Constructor
  EnzoCommBlock() throw()
  { }

  /// Destructor
  ~EnzoCommBlock() throw();

  /// Copy constructor
  EnzoCommBlock(const EnzoCommBlock & EnzoCommBlock) throw();

  /// Assignment operator
  EnzoCommBlock & operator= (const EnzoCommBlock & EnzoCommBlock) throw();

#ifdef CONFIG_USE_CHARM

  /// Initialize a migrated Block
  EnzoCommBlock (CkMigrateMessage *m) 
    : CommBlock (m)
  {
    TRACE("EnzoCommBlock::CkMigrateMessage");
  };

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
  }

#endif
  
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

};

#endif /* ENZO_ENZO_COMM_BLOCK_HPP */

