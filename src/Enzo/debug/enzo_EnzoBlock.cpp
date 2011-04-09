// $Id: enzo_EnzoBlock.cpp 2035 2011-02-28 23:47:31Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoBlock.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:02:02 PST 2011
/// @brief    Implementation of the EnzoBlock class

//======================================================================

#include <stdio.h>

#include "debug.decl.h"

#include "mesh_Block.hpp"
#include "enzo_EnzoBlock.hpp"


EnzoBlock::EnzoBlock(int nx, int ny, int nz,
		     double xm, double ym, double zm,
		     double xp, double yp, double zp,
		     int num_field_blocks) throw()
  : Block(thisIndex.x,thisIndex.y,thisIndex.z,
	  nx,ny,nz,xm,ym,zm,xp,yp,zp,num_field_blocks)
{

  // int i,j;

  // for (i=0; i<MAX_DIMENSION; i++) {
  //   AccelerationField[i] = 0;

  //   for (i=0; i<MAX_DIMENSION; i++) {
  //     AccelerationField[i] = 0;
  //     GridLeftEdge[i] = 0;
  //     GridDimension[i] = 0;
  //     GridStartIndex[i] = 0;
  //     GridEndIndex[i] = 0;
  //     CellWidth[i] = 0;
  //   }

  //   for (j=0; j<MAX_NUMBER_OF_BARYON_FIELDS; j++) {
  //     BaryonField[j] = 0;
  //     OldBaryonField[j] = 0;
  //   }

  // }
  CkPrintf ("EnzoBlock()\n");
  // CANNOT BE INITIALIZED HERE SINCE IT REQUIRES EXTENTS
  //  initialize();
}

