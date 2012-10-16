// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialSedovArray3.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Tue Jan  4 19:30:35 PST 2011
/// @brief    Implementation of an array of Sedov problems, one per Block
///
/// This problem is designed for weak-scaling studies.  Each block
/// contains a small Sedov blast problem.  Problem-specific parameters
/// are used to define the initial values, and whether all blasts go
/// at once or intermitently

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialSedovArray3::EnzoInitialSedovArray3 
(const EnzoConfig * config) throw ()
  : Initial(config->initial_cycle, config->initial_time) 
{
  array_[0] = config->enzo_sedov_array[0];
  array_[1] = config->enzo_sedov_array[1];
  array_[2] = config->enzo_sedov_array[2];
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_CHARM

void EnzoInitialSedovArray3::pup (PUP::er &p)
{
  // NOTE: update whenever attributes change

  TRACEPUP;

  Initial::pup(p);
}
#endif

//----------------------------------------------------------------------

void EnzoInitialSedovArray3::enforce 
(
 Block * block,
 const FieldDescr * field_descr,
 const Hierarchy  * hierarchy
 ) throw()

{

  ASSERT("EnzoInitialSedovArray3",
	 "Block does not exist",
	 block != NULL);

  FieldBlock * field_block = block->field_block();

  ASSERT("EnzoInitialSedovArray3",
	 "Insufficient number of fields",
	 field_descr->field_count() >= 4);

  int index_density      = field_descr->field_id("density");
  int index_velocity_x   = field_descr->field_id("velocity_x");
  int index_velocity_y   = field_descr->field_id("velocity_y");
  int index_velocity_z   = field_descr->field_id("velocity_z");
  int index_total_energy = field_descr->field_id("total_energy");

  enzo_float *  d = (enzo_float *) field_block->field_values(index_density);
  enzo_float * vx = (enzo_float *) field_block->field_values(index_velocity_x);
  enzo_float * vy = (enzo_float *) field_block->field_values(index_velocity_y);
  enzo_float * vz = (enzo_float *) field_block->field_values(index_velocity_z);
  enzo_float * te = (enzo_float *) field_block->field_values(index_total_energy);

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  double xbm,ybm,zbm;
  block->lower(&xbm,&ybm,&zbm);

  double xbp,ybp,zbp;
  block->upper(&xbp,&ybp,&zbp);

  double hx,hy,hz;
  field_block->cell_width(xbm,xbp,&hx,
			  ybm,ybp,&hy,
			  zbm,zbp,&hz);

  // Parameters

  const double sedov_density = 1.0;
  const double sedov_p_in  = 1e-5;
  const double sedov_p_out = 1.0;
  const double sedov_radius = 3.5 * hx;
  const double sedov_radius_2 = sedov_radius*sedov_radius;

  const double sedov_te_in = sedov_p_in  / ((EnzoBlock::Gamma - 1.0) * sedov_density);
  const double sedov_te_out= sedov_p_out / ((EnzoBlock::Gamma - 1.0) * sedov_density);

  ASSERT ("EnzoInitialSedovArray3::enforce",
	  "This problem requires Blocks to be cubical",
	  hx == hy && hy == hz);

  
  int gx,gy,gz;
  field_descr->ghosts(index_density,&gx,&gy,&gz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  // background 

  for (int iz=gz; iz<nz+gz; iz++) {
    double z = zbm + (iz - gz + 0.5)*hz;
    for (int iy=gy; iy<ny+gy; iy++) {
      double y = ybm + (iy - gy + 0.5)*hy;
      for (int ix=gx; ix<nx+gx; ix++) {
	double x = xbm + (ix - gx + 0.5)*hx;
	double r2 = x*x + y*y + z*z;

	int i = INDEX(ix,iy,iz,ngx,ngy);
	d[i]  = sedov_density;
	vx[i] = 0.0;
	vy[i] = 0.0;
	vz[i] = 0.0;
	te[i] = sedov_te_out;

      }
    }
  }

  // array of explosions

  // (kx,ky,kz) index of explosion in domain
  // (x,y,z) position in block
  int kxm = 0;
  int kym = 0;
  int kzm = 0;
  int kxp = array_[0];
  int kyp = array_[1];
  int kzp = array_[2];
  TRACE3 ("SEDOV: %d %d %d",kxp,kyp,kzp);

  double xdm,ydm,zdm;
  hierarchy->lower(&xdm,&ydm,&zdm);
  double xdp,ydp,zdp;
  hierarchy->upper(&xdp,&ydp,&zdp);

  double hxa = (xdp-xdm) / array_[0];
  double hya = (ydp-ydm) / array_[1];
  double hza = (zdp-zdm) / array_[2];

  for (int kz=kzm; kz<kzp; kz++) {
    double zc = hza*(0.5+kz);
    for (int ky=kym; ky<kyp; ky++) {
      double yc = hya*(0.5+ky);
      for (int kx=kxm; kx<kxp; kx++) {
	double xc = hxa*(0.5+kx);
	TRACE3("xc,yc,zc = %f %f %f",xc,yc,zc);

	// (explosion center xc,yc,zc)

	for (int iz=gz; iz<nz+gz; iz++) {
	  double z = zbm + (iz - gz + 0.5)*hz - zc;
	  for (int iy=gy; iy<ny+gy; iy++) {
	    double y = ybm + (iy - gy + 0.5)*hy - yc;
	    for (int ix=gx; ix<nx+gx; ix++) {
	      double x = xbm + (ix - gx + 0.5)*hx - xc;
	      double r2 = x*x + y*y + z*z;

	      int i = INDEX(ix,iy,iz,ngx,ngy);
	      
	      if (r2 < sedov_radius_2) te[i] = sedov_te_in;

	    }
	  }
	}

      }
    }
  }


}
