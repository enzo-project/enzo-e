// $Id: mesh_DataBlock.cpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_DataBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Feb 28 13:22:26 PST 2011
/// @todo     Initialize lower / upper in constructor, removing set_extent)
/// @brief    Implementation of the DataBlock object

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

DataBlock::DataBlock(DataDescr * data_descr,
		     int nx, int ny, int nz,
		     int num_field_blocks) throw ()
  : field_block_()

{ 

  FieldDescr * field_descr = data_descr->field_descr();

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

DataBlock::~DataBlock() throw ()
{ 
  // Deallocate field_block_[]
  for (size_t i=0; i<field_block_.size(); i++) {
    delete field_block_[i];
  }
}

//----------------------------------------------------------------------

DataBlock::DataBlock(const DataBlock & data_block) throw ()
  : field_block_()
/// @param     data_block  Object being copied
{
  copy_(data_block);
}

//----------------------------------------------------------------------

DataBlock & DataBlock::operator = (const DataBlock & data_block) throw ()
/// @param     data_block  Source object of the assignment
/// @return    The target assigned object
{
  copy_(data_block);
  return *this;
}

//----------------------------------------------------------------------

const FieldBlock * DataBlock::field_block (int i) const throw()
{ 
  return field_block_.at(i); 
}

//----------------------------------------------------------------------

FieldBlock * DataBlock::field_block (int i) throw()
{ 
  return field_block_.at(i); 
}

//----------------------------------------------------------------------

void DataBlock::extent
(
 double * lower_x, double * upper_x, 
 double * lower_y, double * upper_y,
 double * lower_z, double * upper_z ) const throw ()
{
  if (lower_x) *lower_x = lower_[0];
  if (lower_y) *lower_y = lower_[1];
  if (lower_z) *lower_z = lower_[2];

  if (upper_x) *upper_x = upper_[0];
  if (upper_y) *upper_y = upper_[1];
  if (upper_z) *upper_z = upper_[2];
}

//----------------------------------------------------------------------

void DataBlock::set_extent
(
 double lower_x, double upper_x,
 double lower_y, double upper_y,
 double lower_z, double upper_z ) throw ()

{
  lower_[0] = lower_x;
  lower_[1] = lower_y;
  lower_[2] = lower_z;

  upper_[0] = upper_x;
  upper_[1] = upper_y;
  upper_[2] = upper_z;
}

//======================================================================

void DataBlock::copy_(const DataBlock & data_block) throw()
{
  UNTESTED("DataBlock::create_");

  // Create a copy of field_block_
  field_block_.resize(data_block.field_block_.size());
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (*(data_block.field_block_[i]));
  }
}

//----------------------------------------------------------------------
