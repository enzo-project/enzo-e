// $Id: method_EnzoFactory.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Method] Declaration of the EnzoFactory class

#include "enzo.hpp"

//----------------------------------------------------------------------

Mesh * EnzoFactory::create_mesh
(
 GroupProcess * group_process,
 int nx,  int ny,  int nz,
 int nbx, int nby, int nbz
 ) throw ()
{
  return new EnzoMesh (this,group_process,nx,ny,nz,nbx,nby,nbz);
}

//----------------------------------------------------------------------

Patch * EnzoFactory::create_patch
(
 Mesh * mesh,
 GroupProcess * group_process,
 int nx,   int ny,  int nz,
 int nbx,  int nby, int nbz
 ) throw()
{
  return new EnzoPatch (mesh,this,group_process,nx,ny,nz,nbx,nby,nbz);
}

//----------------------------------------------------------------------

Block * EnzoFactory::create_block
(
 Patch * patch,
 FieldDescr * field_descr,
 int nx, int ny, int nz,
 int num_field_blocks
 ) throw()
{
  return new EnzoBlock (patch, field_descr, nx,ny,nz,num_field_blocks);
}

