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
 FieldDescr * field_descr,
 int ix, int iy, int iz,
 int nx, int ny, int nz,
 double xm, double ym, double zm,
 double xp, double yp, double zp,
 int num_field_blocks
 ) throw()
{
  return new EnzoBlock (field_descr, 
			ix,iy,iz, 
			nx,ny,nz,
			xm,ym,zm, 
			xp,yp,zp, num_field_blocks);
}

