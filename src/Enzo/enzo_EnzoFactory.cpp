// $Id: method_EnzoFactory.hpp 2009 2011-02-22 19:43:07Z bordner $
// See LICENSE_CELLO file for license and copyright information

/// @file     method_EnzoFactory.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Mar 15 15:29:56 PDT 2011
/// @brief    [\ref Method] Declaration of the EnzoFactory class

#include "enzo.hpp"

//----------------------------------------------------------------------

Block * EnzoFactory::create_block
(
 int ix, int iy, int iz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double hx, double hy, double hz,
 int num_field_blocks
 ) throw()
{
#ifdef CONFIG_USE_CHARM

  ERROR("EnzoFactor::create_block",
	"This function should not be called");
  return 0;

#else

  return new EnzoBlock (ix,iy,iz, 
			nx,ny,nz,
			xm,ym,zm, 
			hx,hy,hz, 
			num_field_blocks);
#endif
}

