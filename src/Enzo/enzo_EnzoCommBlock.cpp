// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoCommBlock
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-11-28
/// @brief    Communication class associated with EnzoBlocks

#include "comm.hpp"

//----------------------------------------------------------------------

EnzoCommBlock::EnzoCommBlock
(
 int ibx, int iby, int ibz,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xpm, double ypm, double zpm, // Patch begin
 double xb, double yb, double zb,    // Block width
 CProxy_Patch proxy_patch,
 int patch_id,
 int patch_rank,
 int num_field_blocks
) throw ()
  :  count_refresh_face_(0)
{
  block_ = new EnzoBlock 
    (ibx,iby,ibz, 
     nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     xb,yb,zb, 
     patch_id,
     patch_rank,
     num_field_blocks,
     this);
}

//----------------------------------------------------------------------

EnzoCommBlock::EnzoCommBlock
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xpm, double ypm, double zpm, // Patch begin
 double xb, double yb, double zb,    // Block width
 CProxy_Patch proxy_patch,
 int patch_id,
 int patch_rank,
 int num_field_blocks) throw ()
{
  block = new EnzoBlock
    (nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     xb,yb,zb, 
     proxy_patch,
     patch_id,
     patch_rank,
     num_field_blocks,
     this);
}

//======================================================================

