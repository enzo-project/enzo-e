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
 Patch        * patch, 
 FieldDescr   * field_descr,
 int nx, int ny, int nz,
 int num_field_blocks) throw ()
  : patch_(patch),
    field_block_()

{ 

  printf ("Block()\n");
  // Initialize field_block_[]
  field_block_.resize(num_field_blocks);
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (field_descr,nx,ny,nz);
  }

  // Initialize extent 

  for (int axis=0; axis<3; axis++) {
    lower_[axis]      = 0.0;
    upper_[axis]      = 1.0;
  }
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

void Block::lower(double * xm, double * ym, double * zm) const throw ()
{
  if (xm) *xm = lower_[0];
  if (ym) *ym = lower_[1];
  if (zm) *zm = lower_[2];
}

//----------------------------------------------------------------------

void Block::set_lower(double xm, double ym, double zm) throw ()
{
  lower_[0] = xm;
  lower_[1] = ym;
  lower_[2] = zm;
}

//----------------------------------------------------------------------

void Block::upper(double * xp, double * yp, double * zp) const throw ()
{
  if (xp) *xp = upper_[0];
  if (yp) *yp = upper_[1];
  if (zp) *zp = upper_[2];
}

//----------------------------------------------------------------------

void Block::set_upper(double xp, double yp, double zp) throw ()
{
  upper_[0] = xp;
  upper_[1] = yp;
  upper_[2] = zp;
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
