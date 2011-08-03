// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#include "mesh.hpp"

//----------------------------------------------------------------------

Mesh * Factory::create_mesh
(
 int nx,  int ny,  int nz,
 int nbx, int nby, int nbz
 ) const throw ()
{
  return new Mesh (this,
		   nx,ny,nz,
		   nbx,nby,nbz);
}

//----------------------------------------------------------------------

Patch * Factory::create_patch
(
 GroupProcess * group_process,
 int nx,   int ny,  int nz,
 int nbx,  int nby, int nbz,
 double xm, double ym, double zm,
 double xp, double yp, double zp
 ) const throw()
{
  return new Patch
    (this,group_process,
     nx,ny,nz,
     nbx,nby,nbz,
     xm,ym,zm,
     xp,yp,zp);

}

//----------------------------------------------------------------------

#ifndef CONFIG_USE_CHARM

Block * Factory::create_block
(
 int ibx, int iby, int ibz,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double xb, double yb, double zb,
 int num_field_blocks
 ) const throw()
{

  return new Block (ibx,iby,ibz, 
		    nbx,nby,nbz,
		    nx,ny,nz,
		    xm,ym,zm, 
		    xb,yb,zb, 
		    num_field_blocks);
}
#endif

//----------------------------------------------------------------------
#ifdef CONFIG_USE_CHARM

CProxy_Block Factory::create_block_array
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double xb, double yb, double zb,
 int num_field_blocks
 ) const throw()
{
  return CProxy_Block::ckNew
    (
     nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     xb,yb,zb, 
     num_field_blocks,
     nbx,nby,nbz);
}
#endif

