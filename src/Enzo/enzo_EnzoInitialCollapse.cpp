// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialCollapse.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    Implementation of an array of Collapse problems, one per Block
///
/// This problem is designed for scaling studies of the Gravity
/// solver.  Each block contains a spherical collapse problem.

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------
#define DEBUG_PERFORMANCE
//----------------------------------------------------------------------

void EnzoInitialCollapse::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,array_,3);
  p | radius_relative_;
  
}

//----------------------------------------------------------------------
void EnzoInitialCollapse::enforce_block
(
 Block * block,
 const FieldDescr * field_descr,
 const ParticleDescr * particle_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  if (!block->is_leaf()) return;

  Timer timer;
  timer.start();

  ASSERT("EnzoInitialCollapse",
	 "Block does not exist",
	 block != NULL);

  Field field = block->data()->field();

  // Get Field parameters
  
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  double bxm,bym,bzm;
  double bxp,byp,bzp;
  double hx,hy,hz;
  int gx,gy,gz;

  // domain extents
  double dxm,dym,dzm;
  double dxp,dyp,dzp;
  hierarchy->lower(&dxm,&dym,&dzm);
  hierarchy->upper(&dxp,&dyp,&dzp);

  // Block extents
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);
  field.cell_width(bxm,bxp,&hx,
		   bym,byp,&hy,
		   bzm,bzp,&hz);
  field.ghost_depth(0,&gx,&gy,&gz);

  const int mx = nx + 2*gx;
  const int my = ny + 2*gy;
  const int mz = nz + 2*gz;

  const int m = mx*my*mz;

  // Get Fields
  
  enzo_float *  d = (enzo_float *) field.values ("density");
  enzo_float * dt = (enzo_float *) field.values ("density_total");
  enzo_float *  p = (enzo_float *) field.values ("pressure");
  enzo_float * po = (enzo_float *) field.values ("potential");
  enzo_float *  t = (enzo_float *) field.values ("temperature");
  enzo_float * te = (enzo_float *) field.values ("total_energy");
  enzo_float * ie = (enzo_float *) field.values ("internal_energy");
  enzo_float * ax = (enzo_float *) field.values ("acceleration_x");
  enzo_float * ay = (enzo_float *) field.values ("acceleration_y");
  enzo_float * az = (enzo_float *) field.values ("acceleration_z");
  enzo_float * vx = (enzo_float *) field.values ("velocity_x");
  enzo_float * vy = (enzo_float *) field.values ("velocity_y");
  enzo_float * vz = (enzo_float *) field.values ("velocity_z");
  enzo_float *  x = (enzo_float *) field.values ("X");
  enzo_float *  b = (enzo_float *) field.values ("B");

  // Initialize Fields

  const int in = cello::index_static();

  const double gamma = EnzoBlock::Gamma[in];
  const double mu = 1.67e-24; // molecular hydrogen
  const double energy = 1e-3*(cello::k)*temperature_ / ((gamma - 1.0) * mu);
  
  // ...compute ellipsoid density

  const double rx = (dxp - dxm) * radius_relative_ / array_[0] ;
  const double ry = (dyp - dym) * radius_relative_ / array_[1] ;
  const double rz = (dzp - dzm) * radius_relative_ / array_[2] ;

  const double rx2i = 1.0/(rx*rx);
  const double ry2i = 1.0/(ry*ry);
  const double rz2i = 1.0/(rz*rz);
  
  const int i0 = gx + mx*(gy + my*gz);
  
  const double density = mass_ / (4.0/3.0*(cello::pi)*rx*ry*rz);

  // bounds of possible explosions intersecting this Block

  int kxm = MAX((int)floor((bxm-dxm-rx)/(dxp-dxm)*array_[0])-1,0);
  int kym = MAX((int)floor((bym-dym-ry)/(dyp-dym)*array_[1])-1,0);
  int kzm = MAX((int)floor((bzm-dzm-rz)/(dzp-dzm)*array_[2])-1,0);
  int kxp = MIN( (int)ceil((bxp-dxm+rx)/(dxp-dxm)*array_[0])+1,array_[0]);
  int kyp = MIN( (int)ceil((byp-dym+ry)/(dyp-dym)*array_[1])+1,array_[1]);
  int kzp = MIN( (int)ceil((bzp-dzm+rz)/(dzp-dzm)*array_[2])+1,array_[2]);
  
  double hxa = (dxp-dxm) / array_[0];
  double hya = (dyp-dym) / array_[1];
  double hza = (dzp-dzm) / array_[2];

  // (kx,ky,kz) index bounds of collapse in domain

  // Initialize background 

  // ratio of density inside and outside the cloud 
  const double density_ratio = 100.0;
  
  for (int i=0; i<m; i++) {
    d[i]  = density / density_ratio;
    te[i] = energy;
    ie[i] = energy;
    t[i]  = temperature_*density_ratio;
    dt[i] = 0.0;
    p[i]  = 0.0;
    po[i] = 0.0;
    ax[i] = 0.0;
    ay[i] = 0.0;
    az[i] = 0.0;
    vx[i] = 0.0;
    vy[i] = 0.0;
    vz[i] = 0.0;
    x[i]  = 0.0;
    b[i]  = 0.0;
  }

  // Initialize sphere (ellipsoid)
  
  for (int kz=kzm; kz<kzp; kz++) {
    double zc = dzm + hza*(0.5+kz);
    for (int ky=kym; ky<kyp; ky++) {
      double yc = dym + hya*(0.5+ky);
      for (int kx=kxm; kx<kxp; kx++) {
	double xc = dxm + hxa*(0.5+kx);
	TRACE3("xc,yc,zc = %f %f %f",xc,yc,zc);

	// (explosion center xc,yc,zc)

	for (int iz=0; iz<mz; iz++) {
	  double z = bzm + (iz - gz + 0.5)*hz - zc;
	  for (int iy=0; iy<my; iy++) {
	    double y = bym + (iy - gy + 0.5)*hy - yc;
	    for (int ix=0; ix<mx; ix++) {
	      double x = bxm + (ix - gx + 0.5)*hx - xc;
	      double r2 = x*x + y*y + z*z;

	      int i = INDEX(ix,iy,iz,mx,my);

	      bool in_sphere = (x*x*rx2i + y*y*ry2i + z*z*rz2i < 1.0);
	      if (in_sphere) {
		d[i]  = density;
                t[i]  = temperature_;
	      }
	    }
	  }
	}
      }
    }
  }

#ifdef DEBUG_PERFORMANCE  
  if (CkMyPe()==0) {
    CkPrintf ("%s:%d %s DEBUG_PERFORMANCE %f\n",
	      __FILE__,__LINE__,block->name().c_str(),
	      timer.value());
  }
#endif  
  // Initialize particles

  Particle particle = block->data()->particle();
  
}

