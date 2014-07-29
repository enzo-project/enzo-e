// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialTurbulence.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:30:49 UTC 2014
/// @brief    Implementation of Enzo 2D Implosion problem initialization

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialTurbulence::EnzoInitialTurbulence 
(int init_cycle, double init_time) throw ()
  : Initial(init_cycle, init_time) 
{ }

//----------------------------------------------------------------------

void EnzoInitialTurbulence::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);
}

//----------------------------------------------------------------------

void EnzoInitialTurbulence::enforce_block 
(
 CommBlock * comm_block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  //  INCOMPLETE("EnzoInitialTurbulence::enforce_block()");

  ASSERT("EnzoInitialTurbulence",
	 "CommBlock does not exist",
	 comm_block != NULL);

  FieldBlock * field_block = comm_block->block()->field_block();

  enzo_float *  d = (enzo_float *) field_block->values("density");
  enzo_float *  p = (enzo_float *) field_block->values("pressure");
  enzo_float * ax = (enzo_float *) field_block->values("driving_x");
  enzo_float * ay = (enzo_float *) field_block->values("driving_y");
  enzo_float * az = (enzo_float *) field_block->values("driving_z");
  enzo_float * vx = (enzo_float *) field_block->values("velocity_x");
  enzo_float * vy = (enzo_float *) field_block->values("velocity_y");
  enzo_float * vz = (enzo_float *) field_block->values("velocity_z");
  enzo_float * te = (enzo_float *) field_block->values("total_energy");

  int rank = comm_block->simulation()->rank();

  ASSERT("EnzoInitializeTurbulence::enforce_block()",
	 "Missing Field 'density'", d);
  ASSERT("EnzoInitializeTurbulence::enforce_block()",
	 "Missing Field 'pressure'",p);
  if (rank >= 1)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'driving_x'", ax);
  if (rank >= 2)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'driving_y'",  ay);
  if (rank >= 3)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'driving_z'",  az);
  if (rank >= 1)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'velocity_x'", vx);
  if (rank >= 2)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'velocity_y'", vy);
  if (rank >= 3)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'velocity_z'",  vz);
  ASSERT("EnzoInitializeTurbulence::enforce_block()",
	 "Missing Field 'total_energy'",  te);

  double xm,ym,zm;
  comm_block->block()->lower(&xm,&ym,&zm);

  double xp,yp,zp;
  comm_block->block()->upper(&xp,&yp,&zp);

  double hx,hy,hz;
  field_block->cell_width(xm,xp,&hx,
			  ym,yp,&hy,
			  zm,zp,&hz);

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);
  int gx,gy,gz;
  field_descr->ghosts(0,&gx,&gy,&gz);

  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;
  int ndz = nz + 2*gz;

  for (int iz=gz; iz<nz+gz; iz++) {
    double z = zm + (iz - gz + 0.5)*hz;
    for (int iy=gy; iy<ny+gy; iy++) {
      double y = ym + (iy - gy + 0.5)*hy;
      for (int ix=gx; ix<nx+gx; ix++) {
	double x = xm + (ix - gx + 0.5)*hx;
	int i = ix + ndx*(iy + ndy*iz);
	d[i]  = 1.0;
	ax[i] = 0.0;
	te[i] = p[i] / ((EnzoBlock::Gamma - 1.0) * d[i]);
      }
    }
  }

  // initialize velocities using turboinit

  int Nx,Ny,Nz;
  comm_block->simulation()->hierarchy()->root_size (&Nx, &Ny, &Nz);

  // assumes cubical domain

  ASSERT3 ("EnzoInitialTurbulence::enforce_block()",
	   "Root grid mesh dimensions %d %d %d must be equal or 1",
	   Nx,Ny,Nz,
	   (Ny==1) || 
	   ((Nz==1) && (Nx == Ny)) ||
	   (Nx == Ny && Ny == Nz));

  int ix,iy,iz;
  comm_block->index().array(&ix,&iy,&iz);

  int ox = ix * nx;
  int oy = iy * ny;
  int oz = iz * nz;

  FORTRAN_NAME(turboinit)
    (&rank, &Nx, (enzo_float *)vx, (enzo_float *)vy, (enzo_float*)vz,
     &ndx,&ndy,&ndz,
     &ox,&oy,&oz);
}
