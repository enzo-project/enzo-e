// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#include "mesh.hpp"

#ifdef CONFIG_USE_CHARM
extern CProxy_Simulation  proxy_simulation;
#endif

//----------------------------------------------------------------------
Hierarchy * Factory::create_hierarchy (int dimension, int refinement) const throw ()
{
  return new Hierarchy (this,dimension,refinement); 
}

//----------------------------------------------------------------------
#ifdef CONFIG_USE_CHARM
CProxy_Patch Factory::create_patch 
#else
Patch * Factory::create_patch 
#endif
(
 const FieldDescr * field_descr,
 int nx,   int ny,  int nz,
 int nx0,  int ny0, int nz0,
 int nbx,  int nby, int nbz,
 double xm, double ym, double zm,
 double xp, double yp, double zp,
 bool allocate_blocks,
 int process_first, int process_last_plus
 ) const throw()
{
#ifdef CONFIG_USE_CHARM
  CProxy_Patch proxy_patch = CProxy_Patch::ckNew
    (nx,ny,nz,
     nx0,ny0,nz0,
     nbx,nby,nbz,
     xm,ym,zm,
     xp,yp,zp,
     allocate_blocks,
     process_first, process_last_plus);
  DEBUG1("proxy_patch.ckLocal() = %p",proxy_patch.ckLocal());
  return proxy_patch;
#else
  return new Patch
    (this,
     field_descr,
     nx,ny,nz,
     nx0,ny0,nz0,
     nbx,nby,nbz,
     xm,ym,zm,
     xp,yp,zp,
     allocate_blocks,
     process_first, process_last_plus);
#endif
}

//----------------------------------------------------------------------

IoBlock * Factory::create_io_block () const throw()
{
  return new IoBlock;
}

//----------------------------------------------------------------------

IoFieldBlock * Factory::create_io_field_block () const throw()
{
  return new IoFieldBlock;

}

//----------------------------------------------------------------------
#ifdef CONFIG_USE_CHARM

CProxy_Block Factory::create_block_array
(
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double xb, double yb, double zb,
 CProxy_Patch proxy_patch,
 int num_field_blocks,
 bool allocate
 ) const throw()
{
  if (allocate) {
    return CProxy_Block::ckNew
      (
       nbx,nby,nbz,
       nx,ny,nz,
       xm,ym,zm, 
       xb,yb,zb, 
       proxy_patch,
       num_field_blocks,
       nbx,nby,nbz);
  } else {
    return CProxy_Block::ckNew();
  }
}

#endif

//----------------------------------------------------------------------
Block * Factory::create_block
(
 int ibx, int iby, int ibz,
 int nbx, int nby, int nbz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double xb, double yb, double zb,
#ifdef CONFIG_USE_CHARM
 CProxy_Patch proxy_patch,
#endif 
 int num_field_blocks
 ) const throw()
{
#ifdef CONFIG_USE_CHARM
    CProxy_Block block_array = CProxy_Block::ckNew
    (nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     xb,yb,zb, 
     proxy_patch,
     num_field_blocks,
     nbx,nby,nbz);
  return block_array(ibx,iby,ibz).ckLocal();
#else
  // CProxy_Block proxy_block_reduce = 
  //   CProxy_Block::ckNew()
  return new Block 
    (ibx,iby,ibz, 
     nbx,nby,nbz,
     nx,ny,nz,
     xm,ym,zm, 
     xb,yb,zb, 
     num_field_blocks);
#endif
}

