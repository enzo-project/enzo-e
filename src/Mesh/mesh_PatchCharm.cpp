// $Id: mesh_PatchCharm.cpp 2181 2011-04-07 00:43:09Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_PatchCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Feb 25 16:20:17 PST 2010
/// @brief    Implementation of the PatchCharm class

#include "cello.hpp"

#include "mesh.hpp"


#ifdef CONFIG_USE_CHARM

//----------------------------------------------------------------------

PatchCharm::PatchCharm
(Factory * factory,
 GroupProcess * group_process,
 int nx,   int ny,  int nz,
 int nbx,  int nby, int nbz,
 double xm, double ym, double zm,
 double xp, double yp, double zp) throw()
  : Patch (factory,group_process,nx,ny,nz,nbx,nby,nbz,xm,ym,zm,xp,yp,zp)
{
}

//----------------------------------------------------------------------

PatchCharm::~PatchCharm() throw()
{
  deallocate_blocks();
}

//----------------------------------------------------------------------

void PatchCharm::allocate_blocks(FieldDescr * field_descr) throw()
{

  // Get number of blocks in the patch

  int nbx,nby,nbz;
  layout_->block_count (&nbx, &nby, &nbz);

  // determine block size

  int mbx = size_[0] / nbx;
  int mby = size_[1] / nby;
  int mbz = size_[2] / nbz;

  // Check that blocks evenly subdivide patch
  if (! ((nbx*mbx == size_[0]) &&
	 (nby*mby == size_[1]) &&
	 (nbz*mbz == size_[2]))) {

    char buffer[ERROR_LENGTH];

    sprintf (buffer,
	     "Blocks must evenly subdivide Patch: "
	     "patch size = (%d %d %d)  block count = (%d %d %d)",
	     size_[0],size_[1],size_[2],
	     nbx,nby,nbz);

    ERROR("PatchCharm::allocate_blocks",  buffer);
      
  }

  // Determine size of each block
  double xb = (upper_[0] - lower_[0]) / nbx;
  double yb = (upper_[1] - lower_[1]) / nby;
  double zb = (upper_[2] - lower_[2]) / nbz;

  CkPrintf ("nbx,nby,nbz = %d %d %d\n",nbx,nby,nbz);

  TRACE("");

  // SEGFAULTS

  block_ = CProxy_BlockCharm::ckNew
    (mbx,mby,mbz,
     lower_[0],lower_[1],lower_[2],
     xb,yb,zb, 1,
     nbx,nby,nbz);

  TRACE("");
}

//----------------------------------------------------------------------

void PatchCharm::deallocate_blocks() throw()
{
  block_.ckDestroy();
}

//----------------------------------------------------------------------

Block * PatchCharm::local_block(size_t i) const throw()
{
  
  ERROR("PatchCharm::local_block",
	"Function should not be called");

  return 0;
}

//======================================================================

#endif /* CONFIG_USE_CHARM */

