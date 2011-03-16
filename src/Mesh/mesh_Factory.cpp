// $Id: mesh_Factory.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Factory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Mesh] Declaration of the Factory class

#include "mesh.hpp"

//----------------------------------------------------------------------

Mesh * Factory::create_mesh
(
 GroupProcess * group_process,
 int nx,  int ny,  int nz,
 int nbx, int nby, int nbz
 ) throw ()
{
  return new Mesh (this,group_process,nx,ny,nz,nbx,nby,nbz);
}

//----------------------------------------------------------------------

Patch * Factory::create_patch
(
 GroupProcess * group_process,
 int nx,   int ny,  int nz,
 int nbx,  int nby, int nbz
 ) throw()
{
  return new Patch (this,group_process,nx,ny,nz,nbx,nby,nbz);
}

//----------------------------------------------------------------------

Block * Factory::create_block
(
 Patch * patch,
 FieldDescr * field_descr,
 int nx, int ny, int nz,
 int num_field_blocks
 ) throw()
{
  return new Block (patch, field_descr, nx,ny,nz,num_field_blocks);
}

