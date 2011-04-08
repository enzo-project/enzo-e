// $Id: mesh_BlockCharm.cpp 2035 2011-02-28 23:47:31Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     mesh_BlockCharm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu Mar  3 23:02:02 PST 2011
/// @brief    Implementation of the BlockCharm class

#include "cello.hpp"

#include "mesh.hpp"

#ifdef CONFIG_USE_CHARM

//======================================================================

BlockCharm::BlockCharm
( int nx, int ny, int nz,
  double xm, double ym, double zm,
  double hx, double hy, double hz,
  int num_field_blocks)
  : Block (thisIndex.x,thisIndex.y,thisIndex.z,
	   nx,ny,nz,
	   xm,ym,zm,
	   hx,hy,hz,num_field_blocks)
{
}


//======================================================================
// CHARM METHODS
//======================================================================

void BlockCharm::p_initial()
{
  CkPrintf ("p_initial %d\n",CkMyPe());
}

// void BlockCharm::p_next()
// {
// }

// void BlockCharm::p_prepare()
// {
// }

// void BlockCharm::p_compute()
// {
// }

#endif /* CONFIG_USE_CHARM */
