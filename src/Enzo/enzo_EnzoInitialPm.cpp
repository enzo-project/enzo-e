// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialPm.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2016-04-29
/// @brief    Implementation of EnzoInitialPm for initializing the PM method

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

void EnzoInitialPm::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  p | field_;
  p | mpp_;
  p | level_;
  p | mask_; // pup::able

}

//----------------------------------------------------------------------

void EnzoInitialPm::enforce_block 
(
 Block * block,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 const Hierarchy  * hierarchy
 ) throw()

{
  Field    field    (block->data()->field());
  Particle particle (block->data()->particle());

  if (mpp_ == 0.0) {
    uniform_placement_ (block,field,particle);
  } else {
    density_placement_ (block,field,particle);
  }
}

//----------------------------------------------------------------------

void EnzoInitialPm::uniform_placement_ 
(Block * block, Field field, Particle particle)
{
  // Insert new tracer particles, one per cell at level_

  // ... get number of cells

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);
  
  const int rank = block->rank();

  // Adjust cell widths for difference in block's level from level_
  const int level_block = block->level();
  double r = 1.0;
  if (level_ != -1) {
    for (int l=level_;      l<level_block; l++) r*=0.5;
    for (int l=level_block; l<level_;      l++) r*=2.0;
  }
  const int rx = (rank >= 1) ? r : 1;
  const int ry = (rank >= 2) ? r : 1;
  const int rz = (rank >= 3) ? r : 1;
  
  // Initialize particle positions

  // ... get block extents

  double xm,ym,zm;
  double xp,yp,zp;
  block->lower(&xm,&ym,&zm);
  block->upper(&xp,&yp,&zp);

  const double xl = xp - xm;
  const double yl = yp - ym;
  const double zl = zp - zm;

  const double hx = xl / nx;
  const double hy = yl / ny;
  const double hz = zl / nz;

  const int it = particle.type_index("dark");

  const int ia_x = (rank >= 1) ? particle.attribute_index (it,"x") : -1;
  const int ia_y = (rank >= 2) ? particle.attribute_index (it,"y") : -1;
  const int ia_z = (rank >= 3) ? particle.attribute_index (it,"z") : -1;

  double * xv = (rank >= 1) ? new double [nx] : NULL;
  double * yv = (rank >= 2) ? new double [ny] : NULL;
  double * zv = (rank >= 3) ? new double [nz] : NULL;

  bool * bitmask = new bool[nx*ny*nz];

  if (rank >= 1) {
    for (int ix=0; ix<nx; ++ix) xv[ix] = xm + (ix + 0.5)*hx;
  }
  if (rank >= 2) {
    for (int iy=0; iy<ny; ++iy) yv[iy] = ym + (iy + 0.5)*hy;
  }
  if (rank >= 3) {
    for (int iz=0; iz<nz; ++iz) zv[iz] = zm + (iz + 0.5)*hz;
  }

  double t = block->time();

  mask_->evaluate (bitmask, t, 
		   nx, nx, xv,
		   ny, ny, yv,
		   nz, nz, zv);

  // count particles
  int np = 0;
  const int rp = (rank <= 1) ? r : (rank <= 2) ? r*r : r*r*r;

  for (int iz=0; iz<nz; ++iz) {
    for (int iy=0; iy<ny; ++iy) {
      for (int ix=0; ix<nx; ++ix) {
	int i = ix + nx*(iy + ny*iz);
	if (bitmask[i]) {
	  np += rp;
	}
      }
    }
  }

  // ... insert uninitialized dark matter particles

  particle.insert_particles (it,np);
  block->simulation()->monitor_insert_particles(np);
  
  const int npb = particle.batch_size();

  int ib=0;  // batch counter
  int ipb=0;  // particle / batch counter 

  double * xa = 0;
  double * ya = 0;
  double * za = 0;

  const int ps  = particle.stride(it,ia_x);

  for (int iz=0; iz<nz; ++iz) {
    for (int iy=0; iy<ny; ++iy) {
      for (int ix=0; ix<nx; ++ix) {
	int i = ix + nx*(iy + ny*iz);
	if (bitmask[i]) {
	  for (int kz = 0; kz<rz; kz++) {
	    for (int ky = 0; ky<ry; ky++) {
	      for (int kx = 0; kx<rx; kx++) {
		if (ipb % npb == 0) {
		  if (rank >= 1) xa = (double *) particle.attribute_array(it,ia_x,ib);
		  if (rank >= 2) ya = (double *) particle.attribute_array(it,ia_y,ib);
		  if (rank >= 3) za = (double *) particle.attribute_array(it,ia_z,ib);
		}
		if (rank >= 1) xa[ipb*ps] = xv[ix] - 0.5*hx + (kx+0.5)*hx/rx;
		if (rank >= 2) ya[ipb*ps] = yv[iy] - 0.5*hy + (ky+0.5)*hy/ry;
		if (rank >= 3) za[ipb*ps] = zv[iz] - 0.5*hz + (kz+0.5)*hz/rz;

		ipb++;

		if (ipb == npb) {
		  ipb=0;
		  ib++;
		}

	      }
	    }
	  }
	}
      }
    }
  }


  delete [] xv;
  delete [] yv;
  delete [] zv;
  delete [] bitmask;

}

//----------------------------------------------------------------------

void EnzoInitialPm::density_placement_ 
(Block * block, Field field, Particle particle)
{

  int did;
  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;

  did = field.field_id( (field_ == "") ? "density" : field_);

  field.dimensions  (did,&mx,&my,&mz);
  field.size           (&nx,&ny,&nz);
  field.ghost_depth (did,&gx,&gy,&gz);

  double * density = (double *) field.values(did);

  // Get cell widths hx,hy,hz

  Data * data = block->data();

  double xm,ym,zm;
  double xp,yp,zp;
  double hx,hy,hz;

  data->lower(&xm,&ym,&zm);
  data->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx, ym,yp,&hy, zm,zp,&hz);

  const int rank = block->rank();

  if (rank < 2) hy = 1.0;
  if (rank < 3) hz = 1.0;

  // ... create vector of integrated mass in the Block

  std::vector<double> ms,xs,ys,zs;
  ms.resize(nx*ny*nz + 1);
  xs.resize(nx*ny*nz + 1);
  ys.resize(nx*ny*nz + 1);
  zs.resize(nx*ny*nz + 1);

  ms[0] = 0.0;
  xs[0] = xm ;
  ys[0] = ym ;
  zs[0] = zm ;
  int ims=1;

  for (int iz=gz; iz<nz+gz; iz++) {
    for (int iy=gy; iy<ny+gy; iy++) {
      for (int ix=gx; ix<nx+gx; ix++) {
	int id = ix + mx*(iy + my*iz);
	double m = density[id] *(hx*hy*hz);
	ms[ims] = ms[ims-1] + m;
	if (rank >= 1) xs[ims-1] = xm + (ix-gx)*hx;
	if (rank >= 2) ys[ims-1] = ym + (iy-gy)*hy;
	if (rank >= 3) zs[ims-1] = zm + (iz-gz)*hz;
	++ ims;
      }
    }
  }

  const double rmax = ms[ims-1];

  // ... compute number of particles to place in the block

  ASSERT1 ("EnzoInitialPm()",
	  "Initial:pm:mpp mass per particle %f must be > 0",
	   mpp_,
	   mpp_ > 0);

  const int np = ms[nx*ny*nz] / mpp_;

  // ... insert uninitialized tracer particles

  const int it = particle.type_index("dark");

  particle.insert_particles (it,np);
  block->simulation()->monitor_insert_particles(np);

  const int ia_x = particle.attribute_index (it,"x");
  const int ia_y = particle.attribute_index (it,"y");
  const int ia_z = particle.attribute_index (it,"z");

  const int npb = particle.batch_size();

  int ib=0;  // batch counter
  int ipb=0;  // particle / batch counter 

  double * xa = 0;
  double * ya = 0;
  double * za = 0;

  const int ps  = particle.stride(it,ia_x);

  for (int ip=0; ip<np; ip++) {

    double r = rmax*rand()/RAND_MAX;

    int imin=0;
    int imax=nx*ny*nz;
    int ims;
    do {
      ims = (imin+imax)/2;
      if (ms[ims] <= r) {
    	imin = ims;
      } else if ( r < ms[ims+1]) {
    	imax = ims;
      }
    } while (imax-imin > 1);
    // for (ims=0; ims<nx*ny*nz-1; ims++) {
    //   if (ms[ims] <= r && r <= ms[ims+1]) break;
    // }
    ims = imin;
    ASSERT6( "EnzoInitialPm",
	     "[%d %d %d] %f <= %f < %f",imin,ims,imax,ms[ims],r,ms[ims+1],
	     (ms[ims] <= r && r < ms[ims+1]));

    //    CkPrintf ("%d %f <= %f < %f\n",ims,ms[ims],r,ms[ims+1]);

    // assert (ims < 0 || ims >= nx*ny*nz) ||
    // ;

    // randomize within cell
    double x = (rank >= 1) ? xs[ims] + hx*rand()/(RAND_MAX+1.0) : 0;
    double y = (rank >= 2) ? ys[ims] + hy*rand()/(RAND_MAX+1.0) : 0;
    double z = (rank >= 3) ? zs[ims] + hz*rand()/(RAND_MAX+1.0) : 0;
    
    // ... if new batch then update position arrays
    if (ipb % npb == 0) {
      if (rank >= 1) xa = (double *) particle.attribute_array (it,ia_x,ib);
      if (rank >= 2) ya = (double *) particle.attribute_array (it,ia_y,ib);
      if (rank >= 3) za = (double *) particle.attribute_array (it,ia_z,ib);
    }

    if (rank >= 1) xa[ipb*ps] = x;
    if (rank >= 2) ya[ipb*ps] = y;
    if (rank >= 3) za[ipb*ps] = z;

    ipb++;

    if (ipb == npb) {
      ipb=0;
      ib++;
    }
   
  }
}
