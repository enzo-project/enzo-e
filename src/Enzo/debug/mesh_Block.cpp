#include <stdio.h>

#include "debug.decl.h"

#include "mesh_Block.hpp"

Block::Block
(
 int ix, int iy, int iz,
 int nx, int ny, int nz,
 double xm, double ym, double zm, // Patch begin
 double xb, double yb, double zb, // Block width
 int num_field_blocks) throw ()

{ 
  // Initialize indices into parent patch

  index_[0] = ix;
  index_[1] = iy;
  index_[2] = iz;

  // Initialize extent 

  lower_[0] = xm + ix*xb;
  lower_[1] = ym + iy*yb;
  lower_[2] = zm + iz*zb;

  upper_[0] = xm + (ix+1)*xb;
  upper_[1] = ym + (iy+1)*yb;
  upper_[2] = zm + (iz+1)*zb;
}
//----------------------------------------------------------------------

Block::~Block()
{ 
}
