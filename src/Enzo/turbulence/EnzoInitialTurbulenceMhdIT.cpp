// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialTurbulenceMhdIT.cpp
/// @author   Alexei Kritsuk (akritsuk@ucsd.edu)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:30:49 UTC 2014
/// @date     Fri Aug 24 00:30:49 UTC 2018
/// @brief    [\ref Enzo] Initial conditions for isothermal MHD turbulence simulation with PPML-IT solver

#include "cello.hpp"

#include "enzo.hpp"

// #define DEBUG_TURBULENCE

#ifdef DEBUG_TURBULENCE
#  define TRACE_TURBULENCE CkPrintf("%s:%d TRACE DEBUG_TURBULENCE\n",__FILE__,__LINE__);
#else
#  define TRACE_TURBULENCE /*  */
#endif
//----------------------------------------------------------------------

EnzoInitialTurbulenceMhdIT::EnzoInitialTurbulenceMhdIT 
(int init_cycle, double init_time, 
 double density_initial,
 double bfieldx_initial,
 double gamma) throw ()
  : Initial(init_cycle, init_time),
    density_initial_(density_initial),
    bfieldx_initial_(bfieldx_initial),
    gamma_(gamma)
{ }

//----------------------------------------------------------------------

void EnzoInitialTurbulenceMhdIT::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  p | density_initial_;
  p | bfieldx_initial_;
  p | gamma_;

}

//----------------------------------------------------------------------

void EnzoInitialTurbulenceMhdIT::enforce_block 
( Block * block, const Hierarchy  * hierarchy ) throw()

{
  TRACE_TURBULENCE;
  if (!block->is_leaf()) return;

  //  INCOMPLETE("EnzoInitialTurbulenceMhdIT::enforce_block()");

  ASSERT("EnzoInitialTurbulenceMhdIT",
	 "Block does not exist",
	 block != NULL);

  Field field = block->data()->field();

  enzo_float * density    = (enzo_float *) field.values("density");
  enzo_float * dens_rx    = (enzo_float *) field.values("dens_rx");
  enzo_float * dens_ry    = (enzo_float *) field.values("dens_ry");
  enzo_float * dens_rz    = (enzo_float *) field.values("dens_rz");

  enzo_float * a3[3]    = { (enzo_float *) field.values("drivx"),
			    (enzo_float *) field.values("drivy"),
			    (enzo_float *) field.values("drivz") };
  enzo_float * a3_rx[3] = { (enzo_float *) field.values("drivx_rx"),
			    (enzo_float *) field.values("drivy_rx"),
			    (enzo_float *) field.values("drivz_rx") };
  enzo_float * a3_ry[3] = { (enzo_float *) field.values("drivx_ry"),
			    (enzo_float *) field.values("drivy_ry"),
			    (enzo_float *) field.values("drivz_ry") };
  enzo_float * a3_rz[3] = { (enzo_float *) field.values("drivx_rz"),
			    (enzo_float *) field.values("drivy_rz"),
			    (enzo_float *) field.values("drivz_rz") };

  enzo_float * v3[3]    = { (enzo_float *) field.values("velox"),
			    (enzo_float *) field.values("veloy"),
			    (enzo_float *) field.values("veloz") };
  enzo_float * v3_rx[3] = { (enzo_float *) field.values("velox_rx"),
			    (enzo_float *) field.values("veloy_rx"),
			    (enzo_float *) field.values("veloz_rx") };
  enzo_float * v3_ry[3] = { (enzo_float *) field.values("velox_ry"),
			    (enzo_float *) field.values("veloy_ry"),
			    (enzo_float *) field.values("veloz_ry") };
  enzo_float * v3_rz[3] = { (enzo_float *) field.values("velox_rz"),
			    (enzo_float *) field.values("veloy_rz"),
			    (enzo_float *) field.values("veloz_rz") };

  enzo_float * b3[3]    = { (enzo_float *) field.values("bfieldx"),
			    (enzo_float *) field.values("bfieldy"),
			    (enzo_float *) field.values("bfieldz") };
  enzo_float * b3_rx[3] = { (enzo_float *) field.values("bfieldx_rx"),
			    (enzo_float *) field.values("bfieldy_rx"),
			    (enzo_float *) field.values("bfieldz_rx") };
  enzo_float * b3_ry[3] = { (enzo_float *) field.values("bfieldx_ry"),
			    (enzo_float *) field.values("bfieldy_ry"),
			    (enzo_float *) field.values("bfieldz_ry") };
  enzo_float * b3_rz[3] = { (enzo_float *) field.values("bfieldx_rz"),
			    (enzo_float *) field.values("bfieldy_rz"),
			    (enzo_float *) field.values("bfieldz_rz") };

  int rank = cello::rank();

  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'density'", density);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'dens_rx'", dens_rx);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'dens_ry'", dens_ry);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'dens_rz'", dens_rz);

  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivx'", a3[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivx_rx'", a3_rx[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivx_ry'", a3_ry[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivx_rz'", a3_rz[0]);
  
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivy'", a3[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivy_rx'", a3_rx[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivy_ry'", a3_ry[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivy_rz'", a3_rz[1]);

  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivz'", a3[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivz_rx'", a3_rx[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivz_ry'", a3_ry[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'drivz_rz'", a3_rz[2]);
  
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'velox'", v3[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'velox_rx'", v3_rx[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'velox_ry'", v3_ry[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'velox_rz'", v3_rz[0]);

  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'veloy'", v3[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'veloy_rx'", v3_rx[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'veloy_ry'", v3_ry[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'veloy_rz'", v3_rz[1]);

  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'veloz'", v3[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'veloz_rx'", v3_rx[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'veloz_ry'", v3_ry[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'veloz_rz'", v3_rz[2]);

  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldx'", b3[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldx_rx'", b3_rx[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldx_ry'", b3_ry[0]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldx_rz'", b3_rz[0]);

  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldy'", b3[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldy_rx'", b3_rx[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldy_ry'", b3_ry[1]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldy_rz'", b3_rz[1]);

  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldz'", b3[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldz_rx'", b3_rx[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldz_ry'", b3_ry[2]);
  ASSERT("EnzoInitializeMHDTurbulenceIT::enforce_block()",
	 "Missing Field 'bfieldz_rz'", b3_rz[2]);


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
  field.ghost_depth(0,&gx,&gy,&gz);

  // ndx(yz) include ghost zones
  // nx(yz) do not

  int ndx = (rank >= 1) ? nx + 2*gx : nx;
  int ndy = (rank >= 2) ? ny + 2*gy : ny;
  int ndz = (rank >= 3) ? nz + 2*gz : nz;

  // initialize driving fields using turboinit

  int Nx,Ny,Nz;
  cello::hierarchy()->root_size (&Nx, &Ny, &Nz);

  // assumes cubical domain

  ASSERT3 ("EnzoInitialTurbulenceMhdIT::enforce_block()",
	   "Root grid mesh dimensions %d %d %d must be equal or 1",
	   Nx,Ny,Nz,
	   ( Ny == 1)  || 
	   ((Nz == 1)  && (Nx == Ny)) ||
	   ((Nx == Ny) && (Ny == Nz)) );

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

  unsigned mask = 1 << (INDEX_BITS_TREE - 1);

  int ix0=gx;
  int iy0=gy;
  int iz0=gz;
  for (int i=0; i<level; i++) {
    ix0 *= 2;
    iy0 *= 2;
    iz0 *= 2;
    bool mx = (mask & bx);
    bool my = (mask & by);
    bool mz = (mask & bz);
    o3[0] = 2*o3[0] + (mx ? nx : 0);
    o3[1] = 2*o3[1] + (my ? ny : 0);
    o3[2] = 2*o3[2] + (mz ? nz : 0);
    mask = mask >> 1;
  }
  o3[0] += ix0 - gx;
  o3[1] += iy0 - gy;
  o3[2] += iz0 - gz;

  // Block zone-centered velocities "velox(yz)" are initialized, including ghost zones

  FORTRAN_NAME(turboinit)
    (&rank, &Nx, 
     (enzo_float *)v3[0],
     (enzo_float *)v3[1],
     (enzo_float *)v3[2],
     &ndx,&ndy,&ndz,
     &o3[0],&o3[1],&o3[2]);

  // set right velocities using linear interpolation

  // Only active block zones initialized for right states
  // Ghost zones for right states will be initialized at refresh step in the preamble of PPML method
  // PPML-IT solver takes 7 primitive variables as input

  for (int i=0; i<ndx*ndy*ndz; i++) {
    for (int id=0; id<rank; id++) {
      v3_rx[id][i] = v3[id][i];
      v3_ry[id][i] = v3[id][i];
      v3_rz[id][i] = v3[id][i];
    }
  }

  /*
  for (int iz=0; iz<nz; iz++) {
    for (int iy=0; iy<ny; iy++) {
      for (int ix=0; ix<nx; ix++) {
	int i = (ix+gx) + ndx*((iy+gy) + ndy*(iz+gz));
	int j = ndx;
	int k = ndx*ndy;
	for (int id=0; id<rank; id++) {
	  v3_rx[id][i] = (v3[id][i]+v3[id][i+1])/2.0;
	  v3_ry[id][i] = (v3[id][i]+v3[id][i+j])/2.0;
	  v3_rz[id][i] = (v3[id][i]+v3[id][i+k])/2.0;
	}
      }
    }
  }
  */
  // set forcing acceleration fields

  for (int i=0; i<ndx*ndy*ndz; i++) {
    for (int id=0; id<rank; id++) {
      a3[id][i] = v3[id][i];
      a3_rx[id][i] = v3_rx[id][i];
      a3_ry[id][i] = v3_ry[id][i];
      a3_rz[id][i] = v3_rz[id][i];
    }
  }

  // set initial densities and B-fields

  for (int i=0; i<ndx*ndy*ndz; i++) {
    density[i]  = density_initial_;
    dens_rx[i]  = density_initial_;
    dens_ry[i]  = density_initial_;
    dens_rz[i]  = density_initial_;
    b3[0][i]    = bfieldx_initial_;
    b3_rx[0][i] = bfieldx_initial_;
    b3_ry[0][i] = bfieldx_initial_;
    b3_rz[0][i] = bfieldx_initial_;
    for (int id=1; id<rank; id++) {
      b3[id][i] = 0.0;
      b3_rx[id][i] = 0.0;
      b3_ry[id][i] = 0.0;
      b3_rz[id][i] = 0.0;
    }
  }

}
