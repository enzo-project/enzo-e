// $Id: enzo_BlockMpi.cpp 2035 2011-02-28 23:47:31Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_BlockMpi.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:02:02 PST 2011
/// @brief    Implementation of the BlockMpi class

#include "cello.hpp"

#include "enzo.hpp"

#ifndef CONFIG_USE_CHARM

//======================================================================

BlockMpi::BlockMpi
( int ix, int iy, int iz,
  int nx, int ny, int nz,
  double xm, double ym, double zm,
  double hx, double hy, double hz,
  int num_field_blocks) throw ()
  : Block (ix,iy,iz,nx,ny,nz,xm,ym,zm,hx,hy,hz,num_field_blocks)
{
}

#endif /* ! CONFIG_USE_CHARM */
