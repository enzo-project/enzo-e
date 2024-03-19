// See LICENSE_CELLO file for license and copyright information

/// @file     problem_MethodTrace.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-11-06
/// @brief    Implementation of the Tracer Particle method

#include "problem.hpp"
#include "charm_simulation.hpp"
  
//----------------------------------------------------------------------

MethodTrace::MethodTrace(ParameterGroup p) noexcept
  : Method ( p.value_float("courant",1.0) ),
    timestep_( p.value_float("timestep", std::numeric_limits<double>::max()) ),
    name_( p.value_string("name", "trace") )
{
  cello::simulation()->refresh_set_name(ir_post_,name_);

  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_particles();
}

//----------------------------------------------------------------------

void MethodTrace::compute ( Block * block) throw()
{
  if (block->is_leaf()) {

    Particle particle (block->data()->particle());
    Field    field    (block->data()->field());

    // initialize trace particle type and position attributes

    const int it = particle.type_index(name_);

    const int ia_x = particle.attribute_index(it,"x");
    const int ia_y = particle.attribute_index(it,"y");
    const int ia_z = particle.attribute_index(it,"z");

    const int dp =  particle.stride(it,ia_x);

    const int rank = cello::rank();

    // get velocity field arrays

    union { float * vxa4; double * vxa8; };
    union { float * vya4; double * vya8; };
    union { float * vza4; double * vza8; };
  
    // NOTE: union so v?a4 also initialized

    vxa8 = (rank >= 1) ? (double *) field.values("velocity_x") : NULL;
    vya8 = (rank >= 2) ? (double *) field.values("velocity_y") : NULL;
    vza8 = (rank >= 3) ? (double *) field.values("velocity_z") : NULL;

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

    double dt = block->state()->dt();

    // Get velocity precision: ASSUMES precision same for all axes and
    // is single or double
    
    const bool is_single =
      (field.precision (field.field_id("velocity_x")) == precision_single);

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
	  double v0 = is_single ? vxa4[ix0] : vxa8[ix0];
	  double v1 = is_single ? vxa4[ix1] : vxa8[ix1];
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

	  const int i00 = ix0+mx*iy0;
	  const int i10 = ix1+mx*iy0;
	  const int i01 = ix0+mx*iy1;
	  const int i11 = ix1+mx*iy1;

	  double vx00 = is_single ? vxa4[i00] : vxa8[i00];
	  double vx10 = is_single ? vxa4[i10] : vxa8[i10];
	  double vx01 = is_single ? vxa4[i01] : vxa8[i01];
	  double vx11 = is_single ? vxa4[i11] : vxa8[i11];

	  double vx = vx00*x1*y1 
	    +         vx10*x0*y1
	    +         vx01*x1*y0
	    +         vx11*x0*y0;

	  double vy00 = is_single ? vya4[i00] : vya8[i00];
	  double vy10 = is_single ? vya4[i10] : vya8[i10];
	  double vy01 = is_single ? vya4[i01] : vya8[i01];
	  double vy11 = is_single ? vya4[i11] : vya8[i11];

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

	  double x1 = 1.0 - x0;
	  double y1 = 1.0 - y0;
	  double z1 = 1.0 - z0;

	  const int i000 = ix0+mx*(iy0 + my*iz0);
	  const int i001 = ix0+mx*(iy0 + my*iz1);
	  const int i010 = ix0+mx*(iy1 + my*iz0);
	  const int i011 = ix0+mx*(iy1 + my*iz1);
	  const int i100 = ix1+mx*(iy0 + my*iz0);
	  const int i101 = ix1+mx*(iy0 + my*iz1);
	  const int i110 = ix1+mx*(iy1 + my*iz0);
	  const int i111 = ix1+mx*(iy1 + my*iz1);

	  double vx000 = is_single ? vxa4[i000] : vxa8[i000];
	  double vx010 = is_single ? vxa4[i010] : vxa8[i010];
	  double vx001 = is_single ? vxa4[i001] : vxa8[i001];
	  double vx011 = is_single ? vxa4[i011] : vxa8[i011];
	  double vx100 = is_single ? vxa4[i100] : vxa8[i100];
	  double vx110 = is_single ? vxa4[i110] : vxa8[i110];
	  double vx101 = is_single ? vxa4[i101] : vxa8[i101];
	  double vx111 = is_single ? vxa4[i111] : vxa8[i111];

	  double vx = vx000*x1*y1*z1
	    +         vx100*x0*y1*z1
	    +         vx010*x1*y0*z1
	    +         vx110*x0*y0*z1
	    +         vx001*x1*y1*z0
	    +         vx101*x0*y1*z0
	    +         vx011*x1*y0*z0
	    +         vx111*x0*y0*z0;

	  double vy000 = is_single ? vya4[i000] : vya8[i000];
	  double vy010 = is_single ? vya4[i010] : vya8[i010];
	  double vy001 = is_single ? vya4[i001] : vya8[i001];
	  double vy011 = is_single ? vya4[i011] : vya8[i011];
	  double vy100 = is_single ? vya4[i100] : vya8[i100];
	  double vy110 = is_single ? vya4[i110] : vya8[i110];
	  double vy101 = is_single ? vya4[i101] : vya8[i101];
	  double vy111 = is_single ? vya4[i111] : vya8[i111];

	  double vy = vy000*x1*y1*z1
	    +         vy100*x0*y1*z1
	    +         vy010*x1*y0*z1
	    +         vy110*x0*y0*z1
	    +         vy001*x1*y1*z0
	    +         vy101*x0*y1*z0
	    +         vy011*x1*y0*z0
	    +         vy111*x0*y0*z0;

	  double vz000 = is_single ? vza4[i000] : vza8[i000];
	  double vz010 = is_single ? vza4[i010] : vza8[i010];
	  double vz001 = is_single ? vza4[i001] : vza8[i001];
	  double vz011 = is_single ? vza4[i011] : vza8[i011];
	  double vz100 = is_single ? vza4[i100] : vza8[i100];
	  double vz110 = is_single ? vza4[i110] : vza8[i110];
	  double vz101 = is_single ? vza4[i101] : vza8[i101];
	  double vz111 = is_single ? vza4[i111] : vza8[i111];

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

//----------------------------------------------------------------------

double MethodTrace::timestep (Block * block) throw()
{
  const int rank = cello::rank();

  double dt = std::numeric_limits<double>::max();

  if (block->is_leaf()) {

    Field    field    = block->data()->field();

    union { float * vxa4; double * vxa8; };
    union { float * vya4; double * vya8; };
    union { float * vza4; double * vza8; };

    // NOTE: union so v?a4 also initialized

    vxa8 = (double *) field.values("velocity_x");
    vya8 = (rank >= 2) ? (double *) field.values("velocity_y") : NULL;
    vza8 = (rank >= 3) ? (double *) field.values("velocity_z") : NULL;

    double xm,ym,zm;
    double xp,yp,zp;
    block->lower(&xm,&ym,&zm);
    block->upper(&xp,&yp,&zp);

    int nx,ny,nz;
    field.size(&nx,&ny,&nz);
    int gx,gy,gz;
    field.ghost_depth(0,&gx,&gy,&gz);
    int mx,my,mz;
    mx = nx + 2*gx;
    my = (rank >= 2) ? ny + 2*gy : 1;
    mz = (rank >= 3) ? nz + 2*gz : 1;

    // Get velocity precision: ASSUMES precision same for all axes and
    // is single or double

    const bool is_single =
      (field.precision (field.field_id("velocity_x")) == precision_single);

    const double hx = (xp-xm)/nx;
    const double hy = (yp-ym)/ny;
    const double hz = (zp-zm)/nz;

    if (rank == 1) {
      for (int ix=0; ix<mx; ix++) {
	int i = ix;
	double vx = is_single ? vxa4[i] : vxa8[i];
	double dt_vx = hx / MAX(fabs(vx),1e-6);
	dt = MIN(dt,dt_vx);
      }
    } else if (rank == 2) {
      for (int iy=0; iy<my; iy++) {
	for (int ix=0; ix<mx; ix++) {
	  int i = ix + mx*iy;
	  double vx = is_single ? vxa4[i] : vxa8[i];
	  double vy = is_single ? vya4[i] : vya8[i];
	  double dt_vx = hx / MAX(fabs(vx),1e-6);
	  double dt_vy = hy / MAX(fabs(vy),1e-6);
	  dt = MIN(dt,dt_vx);
	  dt = MIN(dt,dt_vy);
	}
      }
    } else if (rank == 3) {
      for (int iz=0; iz<mz; iz++) {
	for (int iy=0; iy<my; iy++) {
	  for (int ix=0; ix<mx; ix++) {
	    int i = ix + mx*(iy + my*iz);
	    double vx = fabs(is_single ? vxa4[i] : vxa8[i]);
	    double vy = fabs(is_single ? vya4[i] : vya8[i]);
	    double vz = fabs(is_single ? vza4[i] : vza8[i]);
	    double dt_vx = hx / MAX(fabs(vx),1e-6);
	    double dt_vy = hy / MAX(fabs(vy),1e-6);
	    double dt_vz = hz / MAX(fabs(vz),1e-6);
	    dt = MIN(dt,dt_vx);
	    dt = MIN(dt,dt_vy);
	    dt = MIN(dt,dt_vz);
	  }
	}
      }
    }
  }
  return MIN(timestep_,dt);
}
