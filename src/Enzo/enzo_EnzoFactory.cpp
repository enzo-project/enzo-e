// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Enzo] Declaration of the EnzoFactory class

#include "enzo.hpp"

//----------------------------------------------------------------------
#ifndef CONFIG_USE_CHARM
Block * EnzoFactory::create_block
(
 int ibx, int iby, int ibz,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double hx, double hy, double hz,
 int num_field_blocks
 ) const throw()
{
  return new EnzoBlock 
    (
     ibx,iby,ibz, 
     nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     hx,hy,hz, 
     num_field_blocks);
}
#endif

//----------------------------------------------------------------------
#ifdef CONFIG_USE_CHARM
CProxy_Block EnzoFactory::create_block_array
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double hx, double hy, double hz,
 int num_field_blocks
 ) const throw()
{
  return CProxy_EnzoBlock::ckNew
    (
     nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     hx,hy,hz, 
     num_field_blocks,
     nbx,nby,nbz);
}
#endif


