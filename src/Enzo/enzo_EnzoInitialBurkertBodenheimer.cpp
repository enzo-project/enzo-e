// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialBurkertBodenheimer.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    Implementation of an array of Collapse problems, one per Block
///
/// This problem is designed for scaling studies of the Gravity
/// solver.  Each block contains a spherical collapse problem.

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

// #define DEBUG_PERFORMANCE
#define UNIFORM_DENSITY_PROFILE 1
#define R2_PROFILE              2
//----------------------------------------------------------------------

void EnzoInitialBurkertBodenheimer::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);

  PUParray(p,array_,3);
  p | radius_relative_;

}

//----------------------------------------------------------------------
void EnzoInitialBurkertBodenheimer::enforce_block
( Block * block, const Hierarchy  * hierarchy ) throw()

{

  if (!block->is_leaf()) return;

  Timer timer;
  timer.start();

  ASSERT("EnzoInitialBurkertBodenheimer",
	 "Block does not exist",
	 block != NULL);

  const EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();


  Field field = block->data()->field();

  // Get Field parameters

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  // domain extents
  double dxm,dym,dzm;
  double dxp,dyp,dzp;
  hierarchy->lower(&dxm,&dym,&dzm);
  hierarchy->upper(&dxp,&dyp,&dzp);

  // Block extents
  double bxm,bym,bzm;
  double bxp,byp,bzp;
  block->data()->lower(&bxm,&bym,&bzm);
  block->data()->upper(&bxp,&byp,&bzp);

  const int rank = cello::rank();

  double hx,hy,hz;
  field.cell_width(bxm,bxp,&hx,
		   bym,byp,&hy,
		   bzm,bzp,&hz);
  if (rank < 2) hy = 0.0;
  if (rank < 3) hz = 0.0;


  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  const int mx = nx + 2*gx;
  const int my = ny + 2*gy;
  const int mz = nz + 2*gz;

  const int m = mx*my*mz;

  // Get Fields
  bool RotatingSphere = enzo_config->initial_burkertbodenheimer_rotating;
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

  enzo_float * metal       = field.is_field("metal_density") ?
                             (enzo_float *) field.values("metal_density") : NULL;

  // Initialize Fields

  const int in = cello::index_static();

  const double gamma = EnzoBlock::Gamma[in];
  //const double energy = (1e-3*(cello::kboltz)*temperature_ / ((gamma - 1.0) * (1.0 * cello::mass_hydrogen)))/enzo_units->energy();
  const double energy = (temperature_/enzo_units->temperature()) / ((gamma-1.0)) / enzo_config->ppm_mol_weight;

  // fixed for now about 1 / 10 solar
  const double inner_metal_fraction = 0.0010; // sub-solar
  const double outer_metal_fraction = 1.0E-6; // basically metal free
  // ...compute ellipsoid density

  const double rx = (dxp - dxm) * radius_relative_ / array_[0] ;
  const double ry = (dyp - dym) * radius_relative_ / array_[1] ;
  const double rz = (dzp - dzm) * radius_relative_ / array_[2] ;

  const double rx2i = 1.0/(rx*rx);
  const double ry2i = 1.0/(ry*ry);
  const double rz2i = 1.0/(rz*rz);

  // This is the density at the trucation radius
  const double density = (mass_*cello::mass_solar/enzo_units->mass()) / (4.0/3.0*(cello::pi)*rx*ry*rz);


  // velocity at termination radius (if provided)
  const double outer_velocity = enzo_config->initial_burkertbodenheimer_outer_velocity * 1.0E5 /
                                (enzo_units->length()/enzo_units->time());

/*
  CkPrintf("%s: Density = %e\n", __FUNCTION__, density);
  CkPrintf("%s: mass = %e\n", __FUNCTION__, mass_);
  CkPrintf("%s: rx = %e\n", __FUNCTION__, rx);
  CkPrintf("%s: calculated mass (assuming uniform density) = %e\n",
	 __FUNCTION__, (density*(4.0/3.0*(cello::pi)*rx*ry*rz))* (enzo_units->mass())/(cello::mass_solar));
*/
  //exit(-99);

  // bounds of possible explosions intersecting this Block

  int kxm = MAX((int)floor((bxm-dxm-rx)/(dxp-dxm)*array_[0])-1,0);
  int kym = MAX((int)floor((bym-dym-ry)/(dyp-dym)*array_[1])-1,0);
  int kzm = MAX((int)floor((bzm-dzm-rz)/(dzp-dzm)*array_[2])-1,0);
  int kxp = MIN( (int)ceil((bxp-dxm+rx)/(dxp-dxm)*array_[0])+1,array_[0]);
  int kyp = MIN( (int)ceil((byp-dym+ry)/(dyp-dym)*array_[1])+1,array_[1]);
  int kzp = MIN( (int)ceil((bzp-dzm+rz)/(dzp-dzm)*array_[2])+1,array_[2]);

  double hxa = (dxp-dxm) / array_[0];
  double hya = (rank >= 2) ? (dyp-dym) / array_[1] : 0.0;
  double hza = (rank >= 3) ? (dzp-dzm) / array_[2] : 0.0;

  // (kx,ky,kz) index bounds of collapse in domain

  // Initialize background

  // ratio of density inside and outside the cloud
  const double density_ratio = 228.33;

  std::fill_n(d,m,density / density_ratio);
  std::fill_n(te,m,energy);
  std::fill_n(ie,m,energy);
  std::fill_n(t,m,5000.0);
  std::fill_n(dt,m,0.0);
  std::fill_n(p,m,0.0);
  std::fill_n(po,m,0.0);
  std::fill_n(ax,m,0.0);
  std::fill_n(vx,m,0.0);
  if (metal) std::fill_n(metal,m,outer_metal_fraction*(density/density_ratio));
  if (rank >= 2) std::fill_n(ay,m,0.0);
  if (rank >= 2) std::fill_n(vy,m,0.0);
  if (rank >= 3) std::fill_n(az,m,0.0);
  if (rank >= 3) std::fill_n(vz,m,0.0);
  std::fill_n(x,m,0.0);
  std::fill_n(b,m,0.0);

  // Initialize sphere (ellipsoid)

  for (int kz=kzm; kz<kzp; kz++) {
    double zc = dzm + hza*(0.5+kz);
    for (int ky=kym; ky<kyp; ky++) {
      double yc = dym + hya*(0.5+ky);
      for (int kx=kxm; kx<kxp; kx++) {
	double xc = dxm + hxa*(0.5+kx);

	// (explosion center xc,yc,zc)

	for (int iz=0; iz<mz; iz++) {
	  double z = bzm + (iz - gz + 0.5)*hz - zc;
	  for (int iy=0; iy<my; iy++) {
	    double y = bym + (iy - gy + 0.5)*hy - yc;
	    for (int ix=0; ix<mx; ix++) {
	      double x = bxm + (ix - gx + 0.5)*hx - xc;
	      int i = INDEX(ix,iy,iz,mx,my);
	      double R2 = x*x + y*y + z*z;
              double r2 = x*x*rx2i + y*y*ry2i + z*z*rz2i;
	      bool in_sphere = (r2 < 1.0);
	      if (in_sphere) {
		double radius = sqrt(R2);
		if(R2_PROFILE == densityprofile_) //1/r^2 density profile
		  d[i]  = density*rx*rx/(R2);
		else if(UNIFORM_DENSITY_PROFILE == densityprofile_)
		  d[i] = density;
		else {
		  CkPrintf ("%s:%d %s Unknown densityprofile selected\n",
			    __FILE__,__LINE__,block->name().c_str());
		  CkExit(-99);
		}
		if(RotatingSphere == true) {
		  /* Start with solid body rotation */
		  // Find out which shell the cell is in
      double AngularVelocity = 0.0;


      if (outer_velocity > 0){ // use user-provided value
        AngularVelocity = outer_velocity / (radius_relative_);
      } else { // use prior behavior
		    AngularVelocity =  7.2e-13 * enzo_units->time(); // [rad/s]
      }
		  //float SphereRotationPeriod = 8.72734;
		  //a = Ang(SphereAng1,SphereAng2,rx,sqrt(R2));
		  //double RotVelocityx = -2*(cello::pi)*y / SphereRotationalPeriod;
		  //double RotVelocityy = 2*(cello::pi)*x / SphereRotationalPeriod;
		  //double RotVelocityz = 0.0;
		  vx[i] = -AngularVelocity*y;
		  vy[i] = AngularVelocity*x;
		  vz[i] = 0.0;

		  /* Density perturbation */
		  float cosphi = x/sqrt(x*x+y*y);
		  float sinphi = y/sqrt(x*x+y*y);
		  float phi    = acos(x/sqrt(x*x+y*y));
		  float cos2phi = cosphi*cosphi -sinphi*sinphi;
		  // Burkert & Bodenheimer (1993) m=2 perturbation:
		  float m2mode = 1.0 + 0.1*cos(2.*phi);
		  d[i] += density * m2mode;
		}
                t[i]  = temperature_ / enzo_units->temperature();

    if (metal) metal[i] = d[i]*inner_metal_fraction;

		if(i == 20) {
		  CkPrintf("Density = %e\n", d[i]);
		  CkPrintf("Angular Velocity = %e rad/s", sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i])/radius / enzo_units->time());
		  //CkExit(-99);
		}
	      }
	    }
	  }
	}
      }
    }
  }

//  if (metal){
//    for (int i = 0; i < mx*my*mz; i++) metal[i] = d[i]*outer_metal_fraction;
//  }

  static int counter = 0;
  // CkPrintf("Grand. Block done %d\n", counter++);
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

/************************************************************************/

double Ang(double a1, double a2, double R, double r)
{
  return ((a2-a1)/R)*r + a1;
}
