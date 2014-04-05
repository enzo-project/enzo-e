// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_Block.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-03-10
/// @brief    
///
/// 

#include "mesh.hpp"

//----------------------------------------------------------------------

Block::Block(int nx, int ny, int nz,
	     int num_field_blocks,
	     double xm, double xp,
	     double ym, double yp,
	     double zm, double zp) throw ()
  : num_field_blocks_(num_field_blocks),
    field_block_()
{
  // Initialize field_block_[]
  field_block_.resize(num_field_blocks);
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (nx,ny,nz);
  }
  lower_[0] = xm;
  lower_[1] = ym;
  lower_[2] = zm;
  upper_[0] = xp;
  upper_[1] = yp;
  upper_[2] = zp;
}

//----------------------------------------------------------------------

Block::~Block() throw ()
{
  // Deallocate field_block_[]
  for (size_t i=0; i<field_block_.size(); i++) {
    delete field_block_[i];
    field_block_[i] = 0;
  }
  num_field_blocks_ = 0;
}

//----------------------------------------------------------------------

Block::Block(const Block & block) throw ()
/// @param     block  Object being copied
:
  num_field_blocks_(block.num_field_blocks_)

{
  copy_(block);
}

//----------------------------------------------------------------------

Block & Block::operator= (const Block & block) throw ()
/// @param     block  Source object of the assignment
/// @return    The target assigned object
{
  copy_(block);
  return *this;
}

//----------------------------------------------------------------------

void Block::allocate (const FieldDescr * field_descr) throw()
{
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i]->allocate_array(field_descr,true);
  }
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

void Block::field_cells (double * x, double * y, double * z,
			 int gx, int gy, int gz) const
{
  double xm,ym,zm;
  this->lower(&xm,&ym,&zm);
  double xp,yp,zp;
  this->upper(&xp,&yp,&zp);
  double hx,hy,hz;
  field_block()->cell_width(xm,xp,&hx, ym,yp,&hy, zm,zp,&hz);
    
  int nx,ny,nz;
  field_block()->size(&nx,&ny,&nz);

  int ixm = -gx;
  int iym = -gy;
  int izm = -gz;

  int ixp = nx+gx;
  int iyp = ny+gy;
  int izp = nz+gz;
  
  for (int ix=ixm; ix<ixp; ix++) x[ix-ixm] = xm + (ix+0.5)*hx;
  for (int iy=iym; iy<iyp; iy++) y[iy-iym] = ym + (iy+0.5)*hy;
  for (int iz=izm; iz<izp; iz++) z[iz-izm] = zm + (iz+0.5)*hz;
}

//======================================================================

void Block::copy_(const Block & block) throw()
{
  num_field_blocks_ = block.num_field_blocks_;
  field_block_.resize(block.field_block_.size());
  for (size_t i=0; i<field_block_.size(); i++) {
    field_block_[i] = new FieldBlock (*(block.field_block_[i]));
  }
}
