// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeCicInterp.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-05-05
/// @brief    Implements the EnzoComputeCicInterp class

// #define DEBUG_CIC_INTERP

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputeCicInterp::EnzoComputeCicInterp  
(std::string     field_name,
 std::string     particle_type,
 std::string     particle_attribute,
 double          dt)
  : it_p_ (cello::particle_descr()->type_index (particle_type)),
    ia_p_ (cello::particle_descr()->attribute_index (it_p_,particle_attribute)),
    if_ (cello::field_descr()->field_id (field_name)),
    dt_(dt)
{
}

//----------------------------------------------------------------------

void EnzoComputeCicInterp::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);

  p | it_p_;
  p | ia_p_;
  p | if_;
  p | dt_;
}

//----------------------------------------------------------------------

void EnzoComputeCicInterp::compute ( Block * block) throw()
{

  if (!block->is_leaf()) return;

  compute_(block);
  
}

//----------------------------------------------------------------------

void EnzoComputeCicInterp::compute_(Block * block)
{
  EnzoBlock * enzo_block = enzo::block(block);

  Field field = enzo_block->data()->field();
  Particle particle = enzo_block->data()->particle();

  enzo_float * vf = (enzo_float*)field.values(if_);

  const int ia_x = particle.attribute_position(it_p_,0);
  const int ia_y = particle.attribute_position(it_p_,1);
  const int ia_z = particle.attribute_position(it_p_,2);

  const int ia_vx = particle.attribute_velocity(it_p_,0);
  const int ia_vy = particle.attribute_velocity(it_p_,1);
  const int ia_vz = particle.attribute_velocity(it_p_,2);

  const int dp =  particle.stride(it_p_,ia_x);
  const int da =  particle.stride(it_p_,ia_p_);
  const int dv =  particle.stride(it_p_,ia_vx);

  const int rank = cello::rank();

  int mx,my,mz;
  field.dimensions(0,&mx,&my,&mz);
  const int mxy=mx*my;
  const int i000 = 0;
  const int i001 = mxy;
  const int i010 = mx;
  const int i011 = mx+mxy;
  const int i100 = 1;
  const int i101 = 1+mxy;
  const int i110 = 1+mx;
  const int i111 = 1+mx+mxy;
  
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  // Get block extents and cell widths
  double xm,ym,zm;
  double xp,yp,zp;
  block->lower(&xm,&ym,&zm);
  block->upper(&xp,&yp,&zp);

  const bool lshift = (dt_ != 0.0);
  
  const int nb = particle.num_batches(it_p_);

  if (rank == 1) {

    for (int ib=0; ib<nb; ib++) {

      enzo_float * vp = (enzo_float*) particle.attribute_array(it_p_, ia_p_, ib);

      const int np = particle.num_particles(it_p_,ib);

      enzo_float * xa = (enzo_float *) particle.attribute_array (it_p_,ia_x,ib);
      enzo_float * vxa = lshift ?
	(enzo_float *) particle.attribute_array (it_p_,ia_vx,ib) : nullptr;

      for (int ip=0; ip<np; ip++) {

	enzo_float x = lshift ? xa[ip*dp] + dt_*vxa[ip*dv] : xa[ip*dp];

	enzo_float tx = nx*(x - xm) / (xp - xm) - 0.5;

	int ix0 = gx + floor(tx);

	int ix1 = ix0 + 1;

	enzo_float x0 = 1.0 - (tx - floor(tx));

	enzo_float x1 = 1.0 - x0;

	vp[ip*da] = x0*vf[ix0] + x1*vf[ix1];
      }
    }
  } else if (rank == 2) {

    for (int ib=0; ib<nb; ib++) {

      enzo_float * vp = (enzo_float*) particle.attribute_array(it_p_, ia_p_, ib);

      const int np = particle.num_particles(it_p_,ib);

      enzo_float * xa = (enzo_float *) particle.attribute_array (it_p_,ia_x,ib);
      enzo_float * ya = (enzo_float *) particle.attribute_array (it_p_,ia_y,ib);

      enzo_float * vxa = lshift ?
	(enzo_float *) particle.attribute_array (it_p_,ia_vx,ib) : nullptr;
      enzo_float * vya = lshift ?
	(enzo_float *) particle.attribute_array (it_p_,ia_vy,ib) : nullptr;

      for (int ip=0; ip<np; ip++) {

	enzo_float x = lshift ? xa[ip*dp] + dt_*vxa[ip*dv] : xa[ip*dp];
	enzo_float y = lshift ? ya[ip*dp] + dt_*vya[ip*dv] : ya[ip*dp];

	enzo_float tx = nx*(x - xm) / (xp - xm) - 0.5;
	enzo_float ty = ny*(y - ym) / (yp - ym) - 0.5;

	int ix0 = gx + floor(tx);
	int iy0 = gy + floor(ty);

	int ix1 = ix0 + 1;
	int iy1 = iy0 + 1;

	enzo_float x0 = 1.0 - (tx - floor(tx));
	enzo_float y0 = 1.0 - (ty - floor(ty));

	enzo_float x1 = 1.0 - x0;
	enzo_float y1 = 1.0 - y0;

	enzo_float * vf0 = vf+ix0+mx*iy0;
	
	vp[ip*da] = x0*(y0*vf0[i000] + y1*vf0[i010])
	  +         x1*(y0*vf0[i100] + y1*vf0[i110]);

      }
    }

  } else if (rank == 3) {

    for (int ib=0; ib<nb; ib++) {

      enzo_float * vp = (enzo_float*) particle.attribute_array(it_p_, ia_p_, ib);

      const int np = particle.num_particles(it_p_,ib);

      enzo_float * xa = (enzo_float *) particle.attribute_array (it_p_,ia_x,ib);
      enzo_float * ya = (enzo_float *) particle.attribute_array (it_p_,ia_y,ib);
      enzo_float * za = (enzo_float *) particle.attribute_array (it_p_,ia_z,ib);

      enzo_float * vxa = lshift ?
	(enzo_float *) particle.attribute_array (it_p_,ia_vx,ib) : nullptr;
      enzo_float * vya = lshift ?
	(enzo_float *) particle.attribute_array (it_p_,ia_vy,ib) : nullptr;
      enzo_float * vza = lshift ?
	(enzo_float *) particle.attribute_array (it_p_,ia_vz,ib) : nullptr;
      
      for (int ip=0; ip<np; ip++) {

	enzo_float x = lshift ? xa[ip*dp] + dt_*vxa[ip*dv] : xa[ip*dp];
	enzo_float y = lshift ? ya[ip*dp] + dt_*vya[ip*dv] : ya[ip*dp];
	enzo_float z = lshift ? za[ip*dp] + dt_*vza[ip*dv] : za[ip*dp];

	enzo_float tx = nx*(x - xm) / (xp - xm) - 0.5;
	enzo_float ty = ny*(y - ym) / (yp - ym) - 0.5;
	enzo_float tz = nz*(z - zm) / (zp - zm) - 0.5;

	int ix0 = gx + floor(tx);
	int iy0 = gy + floor(ty);
	int iz0 = gz + floor(tz);

	int ix1 = ix0 + 1;
	int iy1 = iy0 + 1;
	int iz1 = iz0 + 1;

	enzo_float x0 = 1.0 - (tx - floor(tx));
	enzo_float y0 = 1.0 - (ty - floor(ty));
	enzo_float z0 = 1.0 - (tz - floor(tz));

	enzo_float x1 = 1.0 - x0;
	enzo_float y1 = 1.0 - y0;
	enzo_float z1 = 1.0 - z0;

	enzo_float * vf0 = vf + ix0+mx*(iy0+my*iz0);

	vp[ip*da] = x0*(y0*(z0*vf0[i000] + z1*vf0[i001]) +
			y1*(z0*vf0[i010] + z1*vf0[i011])) 
	  +         x1*(y0*(z0*vf0[i100] + z1*vf0[i101]) +
			y1*(z0*vf0[i110] + z1*vf0[i111]));
      }
    }
  }
}  

