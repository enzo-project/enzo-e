// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-11-06
/// @brief    Implementation of the Tracer Particle method

#include "data.hpp"
#include "charm_simulation.hpp"

//----------------------------------------------------------------------

MethodTrace::MethodTrace 
(
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 double courant,
 double timestep
 ) throw() 
  : Method (courant),
    timestep_(timestep)
{
  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);
  refresh(ir)->add_all_particles(particle_descr->num_types());
}

//----------------------------------------------------------------------

void MethodTrace::compute ( Block * block) throw()
{
  if (block->is_leaf()) {

    Particle particle (block->data()->particle());
    Field    field    (block->data()->field());

    // initialize trace particle type and position attributes

    const int it = particle.type_index("trace");

    const int ia_x = particle.attribute_index(it,"x");
    const int ia_y = particle.attribute_index(it,"y");
    const int ia_z = particle.attribute_index(it,"z");

    const int dp =  particle.stride(it,ia_x);

    // get velocity field arrays
    const int rank = block->rank();
    double  * vxa = (rank >= 1) ? 
      (double *)field.values("velocity_x") : NULL;
    double * vya = (rank >= 2) ? 
      (double *)field.values("velocity_y") : NULL;
    double * vza = (rank >= 3) ? 
      (double *)field.values("velocity_z") : NULL;

    int mx,my,mz;
    field.dimensions(0,&mx,&my,&mz);
    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);

    double xm,ym,zm;
    double xp,yp,zp;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);

    const double hx = (xp-xm)/nx;
    const double hy = (yp-ym)/ny;
    const double hz = (zp-zm)/nz;

    double dt = block -> dt();

    // declare particle position arrays

    float * xa = 0;
    float * ya = 0;
    float * za = 0;

    for (int ib=0; ib<particle.num_batches(it); ib++) {
      xa = (float *) particle.attribute_array (it,ia_x,ib);
      ya = (float *) particle.attribute_array (it,ia_y,ib);
      za = (float *) particle.attribute_array (it,ia_z,ib);
      const int np = particle.num_particles(it,ib);

      if (rank == 1) {

	for (int ip=0; ip<np; ip++) {

	  double x = xa[ip*dp];

	  int ix0 = gx + floor((nx-1)*(x - xm) / (xp - xm));

	  int ix1 = ix0 + 1;

	  double x0 = xm + (ix0-gx+0.5)*hx;

	  double x1 = 1 - x0;

	  double v0 = vxa[ix0];

	  double v1 = vxa[ix1];

	  double vx = v0*x1 + v1*x0;

	  xa[ip*dp] += vx*dt;

	}

      } else if (rank == 2) {

	for (int ip=0; ip<np; ip++) {

	  double x = xa[ip*dp];
	  double y = ya[ip*dp];

	  int ix0 = gx + floor((nx-1)*(x - xm) / (xp - xm));
	  int iy0 = gy + floor((ny-1)*(y - ym) / (yp - ym));

	  int ix1 = ix0 + 1;
	  int iy1 = iy0 + 1;

	  double x0 = xm + (ix0-gx+0.5)*hx;
	  double y0 = ym + (iy0-gy+0.5)*hy;

	  double x1 = 1.0 - x0;
	  double y1 = 1.0 - y0;

	  double vx00 = vxa[ix0+mx*iy0];
	  double vx10 = vxa[ix1+mx*iy0];
	  double vx01 = vxa[ix0+mx*iy1];
	  double vx11 = vxa[ix1+mx*iy1];

	  double vx = vx00*x1*y1 
	    +         vx10*x0*y1
	    +         vx01*x1*y0
	    +         vx11*x0*y0;

	  double vy00 = vya[ix0+mx*iy0];
	  double vy10 = vya[ix1+mx*iy0];
	  double vy01 = vya[ix0+mx*iy1];
	  double vy11 = vya[ix1+mx*iy1];

	  double vy = vy00*x1*y1 
	    +         vy10*x0*y1
	    +         vy01*x1*y0
	    +         vy11*x0*y0;

	  xa[ip*dp] += vx*dt;
	  ya[ip*dp] += vy*dt;

	}
      } else if (rank == 3) {
	for (int ip=0; ip<np; ip++) {

	  double x = xa[ip*dp];
	  double y = ya[ip*dp];
	  double z = za[ip*dp];

	  int ix0 = gx + floor((nx-1)*(x - xm) / (xp - xm));
	  int iy0 = gy + floor((ny-1)*(y - ym) / (yp - ym));
	  int iz0 = gz + floor((nz-1)*(z - zm) / (zp - zm));

	  int ix1 = ix0 + 1;
	  int iy1 = iy0 + 1;
	  int iz1 = iz0 + 1;

	  double x0 = xm + (ix0-gx+0.5)*hx;
	  double y0 = ym + (iy0-gy+0.5)*hy;
	  double z0 = zm + (iz0-gz+0.5)*hz;

	  double x1 = 1.0 - xm;
	  double y1 = 1.0 - ym;
	  double z1 = 1.0 - zm;

	  double vx000 = vxa[ix0+mx*(iy0 + my*iz0)];
	  double vx100 = vxa[ix1+mx*(iy0 + my*iz0)];
	  double vx010 = vxa[ix0+mx*(iy1 + my*iz0)];
	  double vx110 = vxa[ix1+mx*(iy1 + my*iz0)];
	  double vx001 = vxa[ix0+mx*(iy0 + my*iz1)];
	  double vx101 = vxa[ix1+mx*(iy0 + my*iz1)];
	  double vx011 = vxa[ix0+mx*(iy1 + my*iz1)];
	  double vx111 = vxa[ix1+mx*(iy1 + my*iz1)];

	  double vx = vx000*x1*y1*z1 
	    +         vx100*x0*y1*z1 
	    +         vx010*x1*y0*z1 
	    +         vx110*x0*y0*z1 
	    +         vx001*x1*y1*z0 
	    +         vx101*x0*y1*z0 
	    +         vx011*x1*y0*z0 
	    +         vx111*x0*y0*z0;

	  double vy000 = vya[ix0+mx*(iy0 + my*iz0)];
	  double vy100 = vya[ix1+mx*(iy0 + my*iz0)];
	  double vy010 = vya[ix0+mx*(iy1 + my*iz0)];
	  double vy110 = vya[ix1+mx*(iy1 + my*iz0)];
	  double vy001 = vya[ix0+mx*(iy0 + my*iz1)];
	  double vy101 = vya[ix1+mx*(iy0 + my*iz1)];
	  double vy011 = vya[ix0+mx*(iy1 + my*iz1)];
	  double vy111 = vya[ix1+mx*(iy1 + my*iz1)];

	  double vy = vy000*x1*y1*z1 
	    +         vy100*x0*y1*z1 
	    +         vy010*x1*y0*z1 
	    +         vy110*x0*y0*z1 
	    +         vy001*x1*y1*z0 
	    +         vy101*x0*y1*z0 
	    +         vy011*x1*y0*z0 
	    +         vy111*x0*y0*z0;

	  double vz000 = vza[ix0+mx*(iy0 + my*iz0)];
	  double vz100 = vza[ix1+mx*(iy0 + my*iz0)];
	  double vz010 = vza[ix0+mx*(iy1 + my*iz0)];
	  double vz110 = vza[ix1+mx*(iy1 + my*iz0)];
	  double vz001 = vza[ix0+mx*(iy0 + my*iz1)];
	  double vz101 = vza[ix1+mx*(iy0 + my*iz1)];
	  double vz011 = vza[ix0+mx*(iy1 + my*iz1)];
	  double vz111 = vza[ix1+mx*(iy1 + my*iz1)];

	  double vz = vz000*x1*y1*z1 
	    +         vz100*x0*y1*z1 
	    +         vz010*x1*y0*z1 
	    +         vz110*x0*y0*z1 
	    +         vz001*x1*y1*z0 
	    +         vz101*x0*y1*z0 
	    +         vz011*x1*y0*z0 
	    +         vz111*x0*y0*z0;

	  xa[ip*dp] += vx*dt;
	  ya[ip*dp] += vy*dt;
	  za[ip*dp] += vz*dt;

	}
      }
    }
  }
  block->compute_done(); 
}
