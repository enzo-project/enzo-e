// See LICENSE_CELLO file for license and copyright information

/// @file     problem_Classname.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-11-06
/// @brief    Implementation of the InitialTrace class

#include "problem.hpp"

//----------------------------------------------------------------------

void InitialTrace::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change
}

//======================================================================

void InitialTrace::enforce_block
 ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
   ) throw()

{
  Field    field    (block->data()->field());
  Particle particle (block->data()->particle());

  // Insert new tracer particles, one per cell

  // ... get number of cells

  int nx,ny,nz;
  field.field_size (0,&nx,&ny,&nz);
  int n = nx*ny*nz;

  // ... insert uninitialized tracer particles

  int it = particle.type_index("trace");
  particle.insert_particles (it,n);

  // Initialize particle positions

  // ... get block extents

  double xm,ym,zm;
  double xp,yp,zp;
  block->lower(&xm,&ym,&zm);
  block->upper(&xp,&yp,&zp);

  int ia_x = particle.attribute_index(it,"x");
  int ia_y = particle.attribute_index(it,"y");
  int ia_z = particle.attribute_index(it,"z");

  const int np = particle.batch_size();

  int ib=0;  // batch counter
  int ip=0;  // particle counter 

  float * xa = 0;
  float * ya = 0;
  float * za = 0;

  const int dp = particle.stride(it,ia_x);

  for (int iz=0; iz<nz; iz++) {
    float z = (nz>1) ? (zm + (zp-zm)*iz/(nz-1)) : 0.0;
    for (int iy=0; iy<ny; iy++) {
      float y = (ny > 1) ? (ym + (yp-ym)*iy/(ny-1)) : 0.0;
      for (int ix=0; ix<nx; ix++) {
	float x = xm + (xp-xm)*ix/(nx-1);

	const int i= ix + nx*(iy + ny*iz);

	// ... if new batch then update position arrays
	if (i % np == 0) {
	  xa = (float *) particle.attribute_array (it,ia_x,ib);
	  ya = (float *) particle.attribute_array (it,ia_y,ib);
	  za = (float *) particle.attribute_array (it,ia_z,ib);
	}

	xa[ip] = x;
	ya[ip] = y;
	za[ip] = z;

	ip+=dp;
	if (ip/dp == np) {
	  ip=0;
	  ib++;
	}
    
      }
    }
  }
}
