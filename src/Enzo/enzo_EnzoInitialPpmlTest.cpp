// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_PpmlTest.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-06-12
/// @brief    Definition of the PpmlTest class

#include "enzo.hpp"

//======================================================================

EnzoInitialPpmlTest::EnzoInitialPpmlTest
(int cycle, double time,
 const EnzoConfig * enzo_config) throw ()
  : Initial(cycle,time)
{
}

//----------------------------------------------------------------------

void EnzoInitialPpmlTest::enforce_block
( Block * block, const Hierarchy * hierarchy ) throw()
{

  // Problem parameters
  // ...center
  double cx=0.5;
  double cy=0.5;
  double cz=0.5;
  // ... radius-squared
  double r2 = 0.05;

  Field field = block->data()->field();

  enzo_float * density    = (enzo_float *) field.values("density");
  enzo_float * velox      = (enzo_float *) field.values("velox");
  enzo_float * veloy      = (enzo_float *) field.values("veloy");
  enzo_float * veloz      = (enzo_float *) field.values("veloz");
  enzo_float * bfieldx    = (enzo_float *) field.values("bfieldx");
  enzo_float * bfieldy    = (enzo_float *) field.values("bfieldy");
  enzo_float * bfieldz    = (enzo_float *) field.values("bfieldz");
  enzo_float * dens_rx    = (enzo_float *) field.values("dens_rx");
  enzo_float * velox_rx   = (enzo_float *) field.values("velox_rx");
  enzo_float * veloy_rx   = (enzo_float *) field.values("veloy_rx");
  enzo_float * veloz_rx   = (enzo_float *) field.values("veloz_rx");
  enzo_float * bfieldx_rx = (enzo_float *) field.values("bfieldx_rx");
  enzo_float * bfieldy_rx = (enzo_float *) field.values("bfieldy_rx");
  enzo_float * bfieldz_rx = (enzo_float *) field.values("bfieldz_rx");
  enzo_float * dens_ry    = (enzo_float *) field.values("dens_ry");
  enzo_float * velox_ry   = (enzo_float *) field.values("velox_ry");
  enzo_float * veloy_ry   = (enzo_float *) field.values("veloy_ry");
  enzo_float * veloz_ry   = (enzo_float *) field.values("veloz_ry");
  enzo_float * bfieldx_ry = (enzo_float *) field.values("bfieldx_ry");
  enzo_float * bfieldy_ry = (enzo_float *) field.values("bfieldy_ry");
  enzo_float * bfieldz_ry = (enzo_float *) field.values("bfieldz_ry");
  enzo_float * dens_rz    = (enzo_float *) field.values("dens_rz");
  enzo_float * velox_rz   = (enzo_float *) field.values("velox_rz");
  enzo_float * veloy_rz   = (enzo_float *) field.values("veloy_rz");
  enzo_float * veloz_rz   = (enzo_float *) field.values("veloz_rz");
  enzo_float * bfieldx_rz = (enzo_float *) field.values("bfieldx_rz");
  enzo_float * bfieldy_rz = (enzo_float *) field.values("bfieldy_rz");
  enzo_float * bfieldz_rz = (enzo_float *) field.values("bfieldz_rz");

  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;

  field.dimensions (0,&mx,&my,&mz);
  field.size         (&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);

  // Block extents
  double xm,ym,zm;
  double xp,yp,zp;
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);

  // Block cell size
  double dx,dy,dz;
  block->cell_width (&dx,&dy,&dz);

  /// NOTE: assumes face densities are on lower faces
  double dx2=0.5*dx;
  double dy2=0.5*dy;
  double dz2=0.5*dz;
  
  for (int iz=0; iz<mz; iz++) {
    double z = zm + (iz - gz + 0.5)*dz;
    double zd  = (z-cz)    *(z-cz);
    double zd2 = (z-cz-dz2)*(z-cz-dz2);
    for (int iy=0; iy<my; iy++) {
      double y = ym + (iy - gy + 0.5)*dy;
      double yd  = (y-cy)    *(y-cy);
      double yd2 = (y-cy-dy2)*(y-cy-dy2);
      for (int ix=0; ix<mx; ix++) {
	double x = xm + (ix - gx + 0.5)*dx;
	double xd  = (x-cx)    *(x-cx);
	double xd2 = (x-cx-dx2)*(x-cx-dx2);

	bool in   = xd + yd + zd  < r2;
	bool in_x = xd2+ yd + zd  < r2;
	bool in_y = xd + yd2+ zd  < r2;
	bool in_z = xd + yd + zd2 < r2;
	
	int i=ix + mx*(iy + my*iz);

	bfieldx[i]    = 10.0;
	bfieldx_rx[i] = 10.0;
	bfieldx_ry[i] = 10.0;
	bfieldx_rz[i] = 10.0;
	bfieldy[i]    = 0.0;
	bfieldy_rx[i] = 0.0;
	bfieldy_ry[i] = 0.0;
	bfieldy_rz[i] = 0.0;
	bfieldz[i]    = 0.0;
	bfieldz_rx[i] = 0.0;
	bfieldz_ry[i] = 0.0;
	bfieldz_rz[i] = 0.0;
	density[i]    = in   ? 1.0 : 0.1;
	dens_rx[i]    = in_x ? 1.0 : 0.1;
	dens_ry[i]    = in_y ? 1.0 : 0.1;
	dens_rz[i]    = in_z ? 1.0 : 0.1;
	velox[i]      = 0.0;
	velox_rx[i]   = 0.0;
	velox_ry[i]   = 0.0;
	velox_rz[i]   = 0.0;
	veloy[i]      = 0.0;
	veloy_rx[i]   = 0.0;
	veloy_ry[i]   = 0.0;
	veloy_rz[i]   = 0.0;
	veloz[i]      = 0.0;
	veloz_rx[i]   = 0.0;
	veloz_ry[i]   = 0.0;
	veloz_rz[i]   = 0.0;
      }
    }
  }
}
