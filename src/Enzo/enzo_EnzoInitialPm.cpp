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
  // Get density field d

  int did;
  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;

  Field    field    = block->data()->field();
  Particle particle = block->data()->particle();

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

  const double xl = xp - xm;
  const double yl = yp - ym;
  const double zl = zp - zm;

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
	xs[ims-1] = xm + (ix-gx)*hx;
	if (rank >= 2) ys[ims-1] = ym + (iy-gy)*hy;
	if (rank >= 3) zs[ims-1] = zm + (iz-gz)*hz;
	++ ims;
      }
    }
  }
  const double rmax = ms[nx*ny*nz];

  // ... compute number of particles to place in the block

  ASSERT1 ("EnzoInitialPm()",
	  "Initial:pm:mpp mass per particle %f must be > 0",
	   mpp_,
	   mpp_ > 0);

  const int np = ms[nx*ny*nz] / mpp_;

  // ... insert uninitialized tracer particles

  const int it = particle.type_index("dark");

  CkPrintf ("Inserting %d particles\n",np);
  particle.insert_particles (it,np);

  const int ia_x = particle.attribute_index (it,"x");
  const int ia_y = particle.attribute_index (it,"y");
  const int ia_z = particle.attribute_index (it,"z");

  const int npb = particle.batch_size();

  int ib=0;  // batch counter
  int ipb=0;  // particle / batch counter 

  double * xa = 0;
  double * ya = 0;
  double * za = 0;

  const int in = CkMyPe() % MAX_NODE_SIZE;

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
    double x = xs[ims] + hx*rand()/(RAND_MAX+1.0);
    double y = (rank >= 2) ? ys[ims] + hy*rand()/(RAND_MAX+1.0) : 0;
    double z = (rank >= 3) ? zs[ims] + hz*rand()/(RAND_MAX+1.0) : 0;
    
    // ... if new batch then update position arrays
    if (ipb % npb == 0) {
      xa = (double *)   particle.attribute_array (it,ia_x,ib);
      if (rank >= 2) ya = (double *)   particle.attribute_array (it,ia_y,ib);
      if (rank >= 3) za = (double *)   particle.attribute_array (it,ia_z,ib);
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
