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
(FieldDescr    * field_descr,
 std::string     field_name,
 ParticleDescr * particle_descr,
 std::string     particle_type,
 std::string     particle_attribute,
 double          dt)
  : it_p_ (particle_descr->type_index (particle_type)),
    ia_p_ (particle_descr->attribute_index (it_p_,particle_attribute)),
    if_ (field_descr->field_id (field_name)),
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
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

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

  const int rank = block->rank();

  int mx,my,mz;
  field.dimensions(0,&mx,&my,&mz);
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

  for (int ib=0; ib<nb; ib++) {

    enzo_float * vp = (enzo_float*) particle.attribute_array(it_p_, ia_p_, ib);

    const int np = particle.num_particles(it_p_,ib);

    if (rank == 1) {

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

    } else if (rank == 2) {

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

	if ( ! ( (0.0 <= x0 && x0 <= 1.0) ||
		 (0.0 <= y0 && y0 <= 1.0) ||
		 (0.0 <= x1 && x1 <= 1.0) ||
		 (0.0 <= y1 && y1 <= 1.0))) {
	  CkPrintf ("ERROR? %s:%d [xy][01] = %f %f  %f %f\n",
		    __FILE__,__LINE__,x0,y0,x1,y1);
	}
	vp[ip*da] = x0*y0*vf[ix0+mx*iy0] 
	  +         x1*y0*vf[ix1+mx*iy0] 
	  +         x0*y1*vf[ix0+mx*iy1] 
	  +         x1*y1*vf[ix1+mx*iy1];

      }

    } else if (rank == 3) {

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

	vp[ip*da] = x0*y0*z0*vf[ix0+mx*(iy0+my*iz0)] 
	  +         x1*y0*z0*vf[ix1+mx*(iy0+my*iz0)] 
	  +         x0*y1*z0*vf[ix0+mx*(iy1+my*iz0)] 
	  +         x1*y1*z0*vf[ix1+mx*(iy1+my*iz0)]
	  +         x0*y0*z1*vf[ix0+mx*(iy0+my*iz1)] 
	  +         x1*y0*z1*vf[ix1+mx*(iy0+my*iz1)] 
	  +         x0*y1*z1*vf[ix0+mx*(iy1+my*iz1)] 
	  +         x1*y1*z1*vf[ix1+mx*(iy1+my*iz1)];

      }
    }
  }
  
}

