// $Id: mesh_Block.cpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Feb 28 13:22:26 PST 2011
/// @todo     Initialize lower / upper in constructor, removing set_extent)
/// @brief    Implementation of the Block object

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

Block::Block
(
 int ix, int iy, int iz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double xp, double yp, double zp,
 int num_field_blocks) throw ()
  : field_block_()

{ 
  // Initialize field_block_[]
  field_block_.resize(num_field_blocks);
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (nx,ny,nz);
  }

  // Initialize indices into parent patch

  index_[0] = ix;
  index_[1] = iy;
  index_[2] = iz;

  // Initialize extent 

  lower_[axis_x] = xm;
  lower_[axis_y] = ym;
  lower_[axis_z] = zm;

  upper_[axis_x] = xp;
  upper_[axis_y] = yp;
  upper_[axis_z] = zp;
}
//----------------------------------------------------------------------

Block::~Block() throw ()
{ 
  // Deallocate field_block_[]
  for (size_t i=0; i<field_block_.size(); i++) {
    delete field_block_[i];
    field_block_[i] = 0;
  }
}

//----------------------------------------------------------------------

Block::Block(const Block & block) throw ()
  : field_block_()
/// @param     block  Object being copied
{
  copy_(block);
}

//----------------------------------------------------------------------

Block & Block::operator = (const Block & block) throw ()
/// @param     block  Source object of the assignment
/// @return    The target assigned object
{
  copy_(block);
  return *this;
}

//----------------------------------------------------------------------

const FieldBlock * Block::field_block (int i) const throw()
{ 
  return field_block_.at(i); 
}

//----------------------------------------------------------------------

FieldBlock * Block::field_block (int i) throw()
{ 
  return field_block_.at(i); 
}

//----------------------------------------------------------------------

void Block::lower(double * x, double * y, double * z) const throw ()
{
  if (x) *x = lower_[0];
  if (y) *y = lower_[1];
  if (z) *z = lower_[2];
}

//----------------------------------------------------------------------

void Block::upper(double * x, double * y, double * z) const throw ()
{
  if (x) *x = upper_[0];
  if (y) *y = upper_[1];
  if (z) *z = upper_[2];
}

//----------------------------------------------------------------------

// void Block::set_lower(double x, double y, double z) throw ()
// {
//   lower_[0] = x;
//   lower_[1] = y;
//   lower_[2] = z;
// }

// //----------------------------------------------------------------------

// void Block::set_upper(double x, double y, double z) throw ()
// {
//   upper_[0] = x;
//   upper_[1] = y;
//   upper_[2] = z;
// }

//----------------------------------------------------------------------

void Block::index_patch (int * ix=0, int * iy=0, int * iz=0) const throw ()
{
  if (ix) (*ix)=index_[0]; 
  if (iy) (*iy)=index_[1]; 
  if (iz) (*iz)=index_[2]; 
}

//----------------------------------------------------------------------

Block * Block::neighbor (axis_enum axis, face_enum face) const throw()
{
}

//----------------------------------------------------------------------

void Block::refresh_ghosts(const FieldDescr * field_descr,
			   int index_field_set) throw()
{
  field_block_[index_field_set]->refresh_ghosts(field_descr);
}

//======================================================================

void Block::copy_(const Block & block) throw()
{
  UNTESTED("Block::create_");

  // Create a copy of field_block_
  field_block_.resize(block.field_block_.size());
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (*(block.field_block_[i]));
  }
}

//----------------------------------------------------------------------
