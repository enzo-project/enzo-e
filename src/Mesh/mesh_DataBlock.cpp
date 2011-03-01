// $Id: mesh_DataBlock.cpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_DataBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Mon Feb 28 13:22:26 PST 2011
/// @brief    Implementation of the DataBlock object

#include "cello.hpp"

#include "mesh.hpp"

//----------------------------------------------------------------------

DataBlock::DataBlock(DataDescr * data_descr,
		     int nx, int ny, int nz,
		     int num_field_blocks) throw ()
  : data_descr_(data_descr),
    field_block_(),
    boundary_face_()

{ 

  FieldDescr * field_descr = data_descr->field_descr();

  // Initialize field_block_[]
  field_block_.resize(num_field_blocks);
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (field_descr,nx,ny,nz);
  }

  // Initialize boundary_face_[][]
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      boundary_face_[axis][face] = false;
    }
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
  : data_descr_(data_block.data_descr_),
    field_block_()
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
void DataBlock::set_boundary_face
(
 bool value,
 face_enum face, 
 axis_enum axis
 ) throw()
{
  // WARNING: recursive
  UNTESTED("DataBlock::set_boundary_face");
  if (face == face_all) {
    set_boundary_face(value,face_lower,axis);
    set_boundary_face(value,face_upper,axis);
  } else if (axis == axis_all) {
    set_boundary_face(value,face,axis_x);
    set_boundary_face(value,face,axis_y);
    set_boundary_face(value,face,axis_z);
  } else {
    boundary_face_[face][axis] = value;
  }
}

//----------------------------------------------------------------------
bool DataBlock::boundary_face(face_enum face,
			       axis_enum axis) throw()
{
  return boundary_face_[face][axis];
}

//======================================================================

void DataBlock::copy_(const DataBlock & data_block) throw()
{
  UNTESTED("DataBlock::create_");

  const FieldBlock * field_block = data_block.field_block();
  const FieldDescr * field_descr = field_block->field_descr();

  // Create a copy of field_block_
  field_block_.resize(data_block.field_block_.size());
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (*field_block_[i]);
  }
  // Copy boundary_face_[][]
  for (int axis=0; axis<3; axis++) {
    for (int face=0; face<2; face++) {
      boundary_face_[axis][face] = data_block.boundary_face_[axis][face];
    }
  }
}

//----------------------------------------------------------------------
