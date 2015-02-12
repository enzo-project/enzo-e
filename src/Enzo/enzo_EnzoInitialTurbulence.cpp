// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialTurbulence.cpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:30:49 UTC 2014
/// @brief    [\ref Enzo] Initial conditions for turbulence simulations

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialTurbulence::EnzoInitialTurbulence 
(int init_cycle, double init_time, 
 double density_initial,
 double pressure_initial,
 double temperature_initial,
 double gamma) throw ()
  : Initial(init_cycle, init_time),
    density_initial_(density_initial),
    pressure_initial_(pressure_initial),
    temperature_initial_(temperature_initial),
    gamma_(gamma)
{ }

//----------------------------------------------------------------------

void EnzoInitialTurbulence::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  p | density_initial_;
  p | pressure_initial_;
  p | temperature_initial_;
  p | gamma_;

}

//----------------------------------------------------------------------

void EnzoInitialTurbulence::enforce_block 
(
 Block * block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  if (!block->is_leaf()) return;

  //  INCOMPLETE("EnzoInitialTurbulence::enforce_block()");

  ASSERT("EnzoInitialTurbulence",
	 "Block does not exist",
	 block != NULL);

  Field field = block->data()->field();

  enzo_float *  d = (enzo_float *) field.values("density");
  enzo_float *  p = (enzo_float *) field.values("pressure");
  enzo_float *  t = (enzo_float *) field.values("temperature");
  enzo_float * a3[3] = { (enzo_float *) field.values("driving_x"),
			 (enzo_float *) field.values("driving_y"),
			 (enzo_float *) field.values("driving_z") };
  enzo_float * v3[3] = { (enzo_float *) field.values("velocity_x"),
			 (enzo_float *) field.values("velocity_y"),
			 (enzo_float *) field.values("velocity_z") };

  enzo_float * te = (enzo_float *) field.values("total_energy");

  int rank = block->simulation()->rank();

  ASSERT("EnzoInitializeTurbulence::enforce_block()",
	 "Missing Field 'density'", d);
  ASSERT("EnzoInitializeTurbulence::enforce_block()",
	 "Missing Field 'pressure'",p);
  ASSERT("EnzoInitializeTurbulence::enforce_block()",
	 "Missing Field 'temperature'",t);
  if (rank >= 1)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'driving_x'", a3[0]);
  if (rank >= 2)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'driving_y'",  a3[1]);
  if (rank >= 3)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'driving_z'",  a3[2]);
  if (rank >= 1)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'velocity_x'", v3[0]);
  if (rank >= 2)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'velocity_y'", v3[1]);
  if (rank >= 3)
    ASSERT("EnzoInitializeTurbulence::enforce_block()",
	   "Missing Field 'velocity_z'",  v3[2]);
  ASSERT("EnzoInitializeTurbulence::enforce_block()",
	 "Missing Field 'total_energy'",  te);

  double xm,ym,zm;
  block->data()->lower(&xm,&ym,&zm);

  double xp,yp,zp;
  block->data()->upper(&xp,&yp,&zp);

  double hx,hy,hz;
  field.cell_width(xm,xp,&hx,
			  ym,yp,&hy,
			  zm,zp,&hz);

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghosts(0,&gx,&gy,&gz);

  int ndx = (rank >= 1) ? nx + 2*gx : nx;
  int ndy = (rank >= 1) ? ny + 2*gy : ny;
  int ndz = (rank >= 1) ? nz + 2*gz : nz;

  // initialize driving fields using turboinit

  int Nx,Ny,Nz;
  block->simulation()->hierarchy()->root_size (&Nx, &Ny, &Nz);

  // assumes cubical domain

  ASSERT3 ("EnzoInitialTurbulence::enforce_block()",
	   "Root grid mesh dimensions %d %d %d must be equal or 1",
	   Nx,Ny,Nz,
	   (Ny==1) || 
	   ((Nz==1) && (Nx == Ny)) ||
	   (Nx == Ny && Ny == Nz));

  // scale by level
  for (int i=0; i<block->level(); i++) Nx  *= 2;

  // compute offsets
  Index index = block->index();

  int ix,iy,iz;
  index.array(&ix,&iy,&iz);

  int bx,by,bz;
  index.tree(&bx,&by,&bz);

  int o3[3] = { ix * nx, iy * ny, iz * nz };

  int level = index.level();

  unsigned mask = 1 << (INDEX_MAX_TREE_BITS - 1);

  for (int i=0; i<level; i++) {
    bool mx = (mask & bx);
    bool my = (mask & by);
    bool mz = (mask & bz);
    o3[0] = 2*o3[0] + (mx ? nx : 0);
    o3[1] = 2*o3[1] + (my ? ny : 0);
    o3[2] = 2*o3[2] + (mz ? nz : 0);
    mask = mask >> 1;
  }

  if ( rank == 2 ) {

    FORTRAN_NAME(turboinit2d)
      (&rank, &Nx, 
       (enzo_float *)v3[0],
       (enzo_float *)v3[1],
       &ndx,&ndy,
       &o3[0],&o3[1]);

  } else if ( rank == 3) {

    FORTRAN_NAME(turboinit)
      (&rank, &Nx, 
       (enzo_float *)v3[0],
       (enzo_float *)v3[1],
       (enzo_float *)v3[2],
       &ndx,&ndy,&ndz,
       &o3[0],&o3[1],&o3[2]);
  }

  if (rank >= 1) for (int i=0; i<ndx*ndy*ndz; i++) a3[0][i] = v3[0][i];
  if (rank >= 2) for (int i=0; i<ndx*ndy*ndz; i++) a3[1][i] = v3[1][i];
  if (rank >= 3) for (int i=0; i<ndx*ndy*ndz; i++) a3[2][i] = v3[2][i];


  for (int iz=0; iz<ndz; iz++) {
    for (int iy=0; iy<ndy; iy++) {
      for (int ix=0; ix<ndx; ix++) {
	int i = ix + ndx*(iy + ndy*iz);
	d[i]  = density_initial_;
      }
    }
  }

  const bool pressure_defined = (pressure_initial_ != 0.0);
  const bool temperature_defined = (temperature_initial_ != 0.0);

  if (pressure_defined) {
    for (int iz=0; iz<ndz; iz++) {
      for (int iy=0; iy<ndy; iy++) {
	for (int ix=0; ix<ndx; ix++) {
	  int i = ix + ndx*(iy + ndy*iz);
	  p[i] = pressure_initial_;
	}
      }
    }
  } else {
    bool comoving_coordinates = false;
    EnzoComputePressure compute_pressure(gamma_,comoving_coordinates);
    compute_pressure.compute(block);
  }

  if (temperature_defined) {
    for (int iz=0; iz<ndz; iz++) {
      for (int iy=0; iy<ndy; iy++) {
	for (int ix=0; ix<ndx; ix++) {
	  int i = ix + ndx*(iy + ndy*iz);
	  te[i] = temperature_initial_;
	  for (int id=0; id<rank; id++) {
	    te[i] += 0.5*v3[id][i]*v3[id][i];
	  }
	}
      }
    }
  } else {
    ERROR("EnzoInitialTurbulence",
	  "Temperature computation with undefined pressure "
	  "not implemented yet");
  }

}
