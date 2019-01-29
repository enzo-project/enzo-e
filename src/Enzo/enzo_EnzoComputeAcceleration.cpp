// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeAcceleration.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoComputeAcceleration class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputeAcceleration::EnzoComputeAcceleration 
(int         rank,
 int         order)
  : Compute(),
    rank_(rank),
    order_(order)
{
  FieldDescr * field_descr = cello::field_descr();
  
  i_ax_ = (rank_ >= 1) ? field_descr->field_id("acceleration_x") : -1;
  i_ay_ = (rank_ >= 2) ? field_descr->field_id("acceleration_y") : -1;
  i_az_ = (rank_ >= 3) ? field_descr->field_id("acceleration_z") : -1;

  i_p_ =  field_descr->field_id("potential");

  if (order_ != 2 && order_ != 4) {
    ERROR1("EnzoComputeAcceleration",
	   "Unknown order %d", order_);
  }
}

//----------------------------------------------------------------------

void EnzoComputeAcceleration::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);

  p | rank_;
  p | order_;
  p | i_ax_;
  p | i_ay_;
  p | i_az_;
  p | i_p_;

}

//----------------------------------------------------------------------

void EnzoComputeAcceleration::compute ( Block * block) throw()
{
  compute_(block);
}

//----------------------------------------------------------------------

void EnzoComputeAcceleration::compute_(Block * block)
{
  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = enzo::block(block);

  Field field = enzo_block->data()->field();

  enzo_float * ax = (rank_ >= 1) ? (enzo_float*) field.values(i_ax_) : NULL;
  enzo_float * ay = (rank_ >= 2) ? (enzo_float*) field.values(i_ay_) : NULL;
  enzo_float * az = (rank_ >= 3) ? (enzo_float*) field.values(i_az_) : NULL;
  enzo_float * p  = (enzo_float*) field.values(i_p_);

  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions  (0,&mx,&my,&mz);
  field.ghost_depth (0,&gx,&gy,&gz);

  int dx,dy,dz;
  dx = 1; 
  dy = mx;
  dz = mx*my;

  int dx2,dy2,dz2;
  dx2 = 2*dx;
  dy2 = 2*dy;
  dz2 = 2*dz;

  double xm,ym,zm;
  double xp,yp,zp;
  double hx,hy,hz;

  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,
		   ym,yp,&hy,
		   zm,zp,&hz);

  // Update cell widths (hx,hy,hz) if needed for expansion

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  enzo_float cosmo_a = 1.0;
  enzo_float cosmo_dadt = 0.0;
  double time = block->time();
  double dt   = block->dt();
  
  if (cosmology) {
   
    cosmology-> compute_expansion_factor (&cosmo_a,&cosmo_dadt,time+0.5*dt);

    hx *= cosmo_a;
    hy *= cosmo_a;
    hz *= cosmo_a;
    dt /= cosmo_a;
  }

  if (order_ == 2) {

    if (rank_ == 1) {

      const enzo_float fx = 1.0 / (2.0*hx);
      for (int ix=1; ix<mx-1; ix++) {
	int i=ix;
	ax[i] = fx*(p[i+dx] - p[i-dx]);
      }

    } else if (rank_ == 2) {

      const enzo_float fx = 1.0 / (2.0*hx);
      const enzo_float fy = 1.0 / (2.0*hy);
      for (int iy=1; iy<my-1; iy++) {
	for (int ix=1; ix<mx-1; ix++) {
	  int i=ix + mx*iy;
	  ax[i] = fx*(p[i+dx] - p[i-dx]);
	  ay[i] = fy*(p[i+dy] - p[i-dy]);
	}
      }

    } else if (rank_ == 3) {

      const enzo_float fx = 1.0 / (2.0*hx);
      const enzo_float fy = 1.0 / (2.0*hy);
      const enzo_float fz = 1.0 / (2.0*hz);
      for (int iz=1; iz<mz-1; iz++) {
	for (int iy=1; iy<my-1; iy++) {
	  for (int ix=1; ix<mx-1; ix++) {
	    int i=ix + mx*(iy + my*iz);
	    ax[i] = fx*(p[i+dx] - p[i-dx]);
	    ay[i] = fy*(p[i+dy] - p[i-dy]);
	    az[i] = fz*(p[i+dz] - p[i-dz]);
	  }
	}
      }
    }

  } else if (order_ == 4) {

    if (rank_ == 1) {

      const enzo_float fx = 1.0 / (12.0*hx);

      for (int ix=2; ix<mx-2; ix++) {
	int i=ix;
	ax[i] = fx*( -p[i+dx2] + 8*p[i+dx] - 8*p[i-dx] + p[i-dx2]);
      }

    } else if (rank_ == 2) {

      const enzo_float fx = 1.0 / (12.0*hx);
      const enzo_float fy = 1.0 / (12.0*hy);
      for (int iy=2; iy<my-2; iy++) {
	for (int ix=2; ix<mx-2; ix++) {
	  int i=ix + mx*iy;
	  ax[i] = fx*( -p[i+dx2] + 8*p[i+dx] - 8*p[i-dx] + p[i-dx2]);
	  ay[i] = fy*( -p[i+dy2] + 8*p[i+dy] - 8*p[i-dy] + p[i-dy2]);
	}
      }

    } else if (rank_ == 3) {

      const enzo_float fx = 1.0 / (12.0*hx);
      const enzo_float fy = 1.0 / (12.0*hy);
      const enzo_float fz = 1.0 / (12.0*hz);
      for (int iz=2; iz<mz-2; iz++) {
	for (int iy=2; iy<my-2; iy++) {
	  for (int ix=2; ix<mx-2; ix++) {
	    int i=ix + mx*(iy + my*iz);
	    ax[i] = fx*( -p[i+dx2] + 8*p[i+dx] - 8*p[i-dx] + p[i-dx2]);
	    ay[i] = fy*( -p[i+dy2] + 8*p[i+dy] - 8*p[i-dy] + p[i-dy2]);
	    az[i] = fz*( -p[i+dz2] + 8*p[i+dz] - 8*p[i-dz] + p[i-dz2]);
	  }
	}
      }

    }

  } else {
    ERROR1("EnzoComputeAcceleration",
	   "Unknown order %d", order_);
  }

  // Update particle accelerations

  Particle particle = block->data()->particle();

  ParticleDescr * particle_descr = cello::particle_descr();
  Grouping * particle_groups = particle_descr->groups();

  int num_mass = particle_groups->size("has_mass");

  for (int ipt = 0; ipt < num_mass; ipt++){

    std::string particle_type = particle_groups->item("has_mass",ipt);
    int it = particle.type_index(particle_type);

    if (particle.num_particles(it) > 0) {

      double dt_shift = 0.5*block->dt() / cosmo_a;
      //  double dt_shift = 0.0;
      if (rank_ >= 1) {
        EnzoComputeCicInterp interp_x
        	("acceleration_x", particle_type, "ax",dt_shift);
        interp_x.compute(block);
      }
      if (rank_ >= 2) {
        EnzoComputeCicInterp interp_y
	        ("acceleration_y", particle_type, "ay",dt_shift);
        interp_y.compute(block);
      }
      if (rank_ >= 3) {
        EnzoComputeCicInterp interp_z
	        ("acceleration_z", particle_type, "az",dt_shift);
        interp_z.compute(block);
      }
    }
  }
}

