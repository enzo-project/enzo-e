// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulenceMhdIG.cpp
/// @author   Alexei Kritsuk (akritsuk@ucsd.edu)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:31:04 UTC 2014
/// @date     Thu Sep 20 00:31:04 UTC 2018
/// @brief    Implements the EnzoMethodTurbulenceMhdIG class with Ornstein-Uhlenbeck pumping

#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

enum {
  id_vad,
  id_aad,
  id_vvd,
  id_vv,
  id_dvx,
  id_dvy,
  id_dvz,
  id_dax,
  id_day,
  id_daz,
  id_bx,
  id_by,
  id_bz,
  id_bb,
  id_bbod,
  id_divb,
  id_d,
  id_dd,
  id_lnd,
  id_dlnd,
  id_pr,
  id_prod,
  id_zones,
  id_mind,
  id_maxd,
  num_reduce
};

//#define DEBUG_TURBULENCE

#ifdef DEBUG_TURBULENCE
#   define TRACE_TURBULENCE CkPrintf ("%s:%d TRACE DEBUG_TURBULENCE\n",__FILE__,__LINE__);
#else
#   define TRACE_TURBULENCE /*   */
#endif

//----------------------------------------------------------------------

CkReduction::reducerType r_method_turbulence_id_type;

void register_method_turbulence_mhd_ig(void)
{ 
  r_method_turbulence_id_type = CkReduction::addReducer(r_method_turbulence_mhd_ig); 
}

CkReductionMsg * r_method_turbulence_mhd_ig(int n, CkReductionMsg ** msgs)
{
  double accum[num_reduce] = { 0.0 };
  accum[id_mind] =   std::numeric_limits<double>::max();
  accum[id_maxd] = - std::numeric_limits<double>::max();

  for (int i=0; i<n; i++) {
    double * values = (double *) msgs[i]->getData();
    for (int ig=0; ig<num_reduce-2; ig++) {
      accum [ig] += values[ig];
    }
    accum [id_mind] = 
      std::min(accum[id_mind],values[id_mind]);
    accum [id_maxd] = 
      std::max(accum[id_maxd],values[id_maxd]);
  }
  return CkReductionMsg::buildNew(num_reduce*sizeof(double),accum);
}

//----------------------------------------------------------------------

EnzoMethodTurbulenceMhdIG::EnzoMethodTurbulenceMhdIG 
(double gamma,
 double density_initial,
 double pressure_initial,
 double bfieldx_initial,
 double mach_number,
 double solenoidal_fraction,
 double kfmin,
 double kfmax,
 bool comoving_coordinates)
  : Method(),
    gamma_(gamma),
    density_initial_(density_initial),
    pressure_initial_(pressure_initial),
    bfieldx_initial_(bfieldx_initial),
    mach_number_(mach_number),
    solenoidal_fraction_(solenoidal_fraction),
    kfmin_(kfmin),
    kfmax_(kfmax),
    comoving_coordinates_(comoving_coordinates)
{
  TRACE_TURBULENCE;  
  
  // Initialize default Refresh object

  Refresh * refresh_post = cello::refresh(ir_post_);
  cello::simulation()->refresh_set_name(ir_post_,name());
  refresh_post->add_all_fields();

   // TURBULENCE parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceMhdIG::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | gamma_;
  p | density_initial_;
  p | pressure_initial_;
  p | bfieldx_initial_;
  p | mach_number_;
  p | solenoidal_fraction_;
  p | kfmin_;
  p | kfmax_;
  p | comoving_coordinates_;

}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceMhdIG::compute ( Block * block) throw()
{
  TRACE_TURBULENCE;  

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  
  Field field = block->data()->field();

  enzo_float *  density = (enzo_float *) field.values("density");
  enzo_float *  dens_rx = (enzo_float *) field.values("dens_rx");
  enzo_float *  dens_ry = (enzo_float *) field.values("dens_ry");
  enzo_float *  dens_rz = (enzo_float *) field.values("dens_rz");

  enzo_float *  pressure = (enzo_float *) field.values("pressure");
  enzo_float *  press_rx = (enzo_float *) field.values("press_rx");
  enzo_float *  press_ry = (enzo_float *) field.values("press_ry");
  enzo_float *  press_rz = (enzo_float *) field.values("press_rz");

  enzo_float *  velocity[3] = {
    (enzo_float *) field.values("velox"),
    (enzo_float *) field.values("veloy"),
    (enzo_float *) field.values("veloz") };

  enzo_float *  velo_rx[3] = {
    (enzo_float *) field.values("velox_rx"),
    (enzo_float *) field.values("veloy_rx"),
    (enzo_float *) field.values("veloz_rx") };

  enzo_float *  velo_ry[3] = {
    (enzo_float *) field.values("velox_ry"),
    (enzo_float *) field.values("veloy_ry"),
    (enzo_float *) field.values("veloz_ry") };

  enzo_float *  velo_rz[3] = {
    (enzo_float *) field.values("velox_rz"),
    (enzo_float *) field.values("veloy_rz"),
    (enzo_float *) field.values("veloz_rz") };

  enzo_float * driving[3] = {
    (enzo_float *) field.values("drivx"),
    (enzo_float *) field.values("drivy"),
    (enzo_float *) field.values("drivz") };

  enzo_float * driv_rx[3] = {
    (enzo_float *) field.values("drivx_rx"),
    (enzo_float *) field.values("drivy_rx"),
    (enzo_float *) field.values("drivz_rx") };

  enzo_float * driv_ry[3] = {
    (enzo_float *) field.values("drivx_ry"),
    (enzo_float *) field.values("drivy_ry"),
    (enzo_float *) field.values("drivz_ry") };

  enzo_float * driv_rz[3] = {
    (enzo_float *) field.values("drivx_rz"),
    (enzo_float *) field.values("drivy_rz"),
    (enzo_float *) field.values("drivz_rz") };

  enzo_float * bfield[3] = {
    (enzo_float *) field.values("bfieldx"),
    (enzo_float *) field.values("bfieldy"),
    (enzo_float *) field.values("bfieldz") };

  enzo_float * bfield_rx[3] = {
    (enzo_float *) field.values("bfieldx_rx"),
    (enzo_float *) field.values("bfieldy_rx"),
    (enzo_float *) field.values("bfieldz_rx") };

  enzo_float * bfield_ry[3] = {
    (enzo_float *) field.values("bfieldx_ry"),
    (enzo_float *) field.values("bfieldy_ry"),
    (enzo_float *) field.values("bfieldz_ry") };

  enzo_float * bfield_rz[3] = {
    (enzo_float *) field.values("bfieldx_rz"),
    (enzo_float *) field.values("bfieldy_rz"),
    (enzo_float *) field.values("bfieldz_rz") };

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);

  const int n = num_reduce;
  double * g = new double [n+1];
  
  for (int i=0; i<n-2; i++) g[i] = 0.0;

  g[id_mind] =   std::numeric_limits<double>::max();
  g[id_maxd] = - std::numeric_limits<double>::max();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  //---------------------------------------------------------------------
  // OU pumping call starts here

  double xdm,ydm,zdm;
  cello::hierarchy()->lower(&xdm,&ydm,&zdm);
  double xdp,ydp,zdp; 
  cello::hierarchy()->upper(&xdp,&ydp,&zdp);

  double Lbox[3];
  Lbox[0] =               (xdp - xdm);
  Lbox[1] = (rank >= 2) ? (ydp - ydm) : 1.0;
  Lbox[2] = (rank >= 3) ? (zdp - zdm) : 1.0;


  // if ( rank == 3) {

  //   FORTRAN_NAME(OUpumpInit)
  //     ( &gamma_, 
  //   	&density_initial_,
  //   	&pressure_initial_,
  //   	&solenoidal_fraction_, 
  //   	&mach_number_, 
  //   	&kfmin_, 
  //   	&kfmax_,
  //   	Lbox );


  //   // FORTRAN_NAME(OUpumpCompute)
  //   //   ( &rank, &mx,  &my,  &mz,         // rank and local block dimensions
  //   //     &nx,  &ny, &nz,              // root 
  //   //     &gx,   &gy,  &gz,               // number of ghost zones
  //   //     (enzo_float *)v3[0],            // flow fields invilved
  //   //     (enzo_float *)v3[1],
  //   //     (enzo_float *)v3[2],
  //   //     (enzo_float *)density,
  //   //     &mx,&my,&mz,                 // zone sizes
  //   //     &o3[0],&o3[1],&o3[2],
  //   //     &dt );                          // time step

  //   turbForce3D(nc, ni, nj, nk, nig, njg, nkg, w, grid, dt, res, update_sol)

  // }

  // OU pumping call ends here
  //---------------------------------------------------------------------


  if (block->is_leaf()) {

    // Loops cover only active zones of this block and
    // compute averages for zone centers (mostly) to be used for forcing normalization

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {

	  int i = (ix+gx) + mx*((iy+gy) + my*(iz+gz));

	  enzo_float d   = density[i];
	  enzo_float pr  = pressure[i];
	  for (int id=0; id<rank; id++) {
	    enzo_float v  = velocity[id][i];
	    enzo_float v2 = v*v;
	    enzo_float a  = driving[id][i];
	    enzo_float b  = bfield[id][i];
	    enzo_float b2 = b*b;

	    g[id_vad]   += v*a*d; // 0
	    g[id_aad]   += a*a*d; // 1
	    g[id_vvd]   += v2*d;  // 2
	    g[id_vv]    += v2;    // 3
	    g[id_bb]    += b2;    // 13
	    g[id_bbod]  += b2/d;  // 14
	  }
	  g[id_d]    +=  d;       // 16
	  g[id_dd]   +=  d*d;     // 17
	  g[id_lnd]  +=  log(d);  // 18
	  g[id_dlnd] +=  d*log(d);// 19

	  g[id_pr]   +=  pr;      // 20
	  g[id_prod] +=  pr/d;    // 21

          g[id_zones] += 1;       // 22

	  g[id_dvx]  +=  d*velocity[0][i];                     // 4
	  g[id_dvy]  +=  (rank >= 2) ? d*velocity[1][i] : 0.0; // 5
	  g[id_dvz]  +=  (rank >= 3) ? d*velocity[2][i] : 0.0; // 6
	  /*
	  g[id_dax]  +=  d*driving[0][i];                      // 7
	  g[id_day]  +=  (rank >= 2) ? d*driving[1][i] : 0.0;  // 8
	  g[id_daz]  +=  (rank >= 3) ? d*driving[2][i] : 0.0;  // 9
	  */
	  g[id_dax]  +=  ( d*driving[0][i] + dens_rx[i]*driv_rx[0][i] + 
	                                                          dens_ry[i]*driv_ry[0][i] + 
						                  dens_rz[i]*driv_rz[0][i] )/4.0;                      // 7
	  g[id_day]  +=  (rank >= 2) ? ( d*driving[1][i] + dens_rx[i]*driv_rx[1][i] +
                                                                                dens_ry[i]*driv_ry[1][i] +
	                                                                        dens_rz[i]*driv_rz[1][i] )/4.0 : 0.0;  // 8
	  g[id_daz]  +=  (rank >= 3) ? ( d*driving[2][i] + dens_rx[i]*driv_rx[2][i] +
                                                                                dens_ry[i]*driv_ry[2][i] +
	                                                                        dens_rz[i]*driv_rz[2][i] )/4.0 : 0.0;  // 9

	  g[id_bx]   +=  bfield[0][i];                         // 10
	  g[id_by]   +=  (rank >= 2) ? bfield[1][i] : 0.0;     // 11
	  g[id_bz]   +=  (rank >= 3) ? bfield[2][i] : 0.0;     // 12

	  // PPML-style divergence calculation in 3D.
	  // One can also use max|div(b)| for control, see Fig. 18 in Ustyugov et al. (2009, JCP 228, 7614).
	  int is = 1;
	  int js = mx;
	  int ks = mx*my;
	  g[id_divb] +=  fabs(bfield_rx[0][i+is] - bfield_rx[0][i-is] +          // 15
						   bfield_rx[0][i+is+js] - bfield_rx[0][i-is+js] +
						   bfield_rx[0][i+is+ks] - bfield_rx[0][i-is+ks] +
						   bfield_rx[0][i+is+js+ks] - bfield_rx[0][i-is+js+ks] +

	                                           bfield_ry[1][i+js] - bfield_ry[1][i-js] +
						   bfield_ry[1][i+js+is] - bfield_ry[1][i-js+is] +
						   bfield_ry[1][i+js+ks] - bfield_ry[1][i-js+ks] +
						   bfield_ry[1][i+js+is+ks] - bfield_ry[1][i-js+is+ks] +

	                                           bfield_rz[2][i+ks] - bfield_rz[2][i-ks] +
						   bfield_rz[2][i+ks+is] - bfield_rz[2][i-ks+is] +
						   bfield_rz[2][i+ks+js] - bfield_rz[2][i-ks+js] +
						   bfield_rz[2][i+ks+is+js] - bfield_rz[2][i-ks+is+js]);

	  g[id_mind] =                                         // 23
	    std::min(g[id_mind], (double) d);
	  g[id_maxd] =                                         // 24
	    std::max(g[id_maxd], (double) d);
	}
      }
    }
  }
  TRACE_TURBULENCE;  
  CkCallback callback (CkIndex_EnzoBlock::r_method_turbulence_ig_end(NULL),
		       enzo_block->proxy_array());
  enzo_block->contribute(n,g,r_method_turbulence_ig_type,callback);
}

//----------------------------------------------------------------------

// SEE main.cpp for implementation

// CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs)

//----------------------------------------------------------------------

void EnzoBlock::r_method_turbulence_ig_end(CkReductionMsg * msg)
{
  TRACE_TURBULENCE;  
  performance_start_(perf_compute,__FILE__,__LINE__);
  method()->compute_resume (this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceMhdIG::compute_resume 
(Block * block,
 CkReductionMsg * msg) throw()
{
  TRACE_TURBULENCE;  

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  double * g = msg->getData();

  Data * data = block->data();
  Field field = data->field();
  
  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  double dt = block->dt();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  double xdm,ydm,zdm;
  cello::hierarchy()->lower(&xdm,&ydm,&zdm);
  double xdp,ydp,zdp; 
  cello::hierarchy()->upper(&xdp,&ydp,&zdp);

  double bnotx    = bfieldx_initial_;

  if (edot_ < 0.0) {
    // Only compute if needed at the beginning--could/should be in
    // EnzoInitialTurbulence
    double domain_x =               (xdp - xdm);
    double domain_y = (rank >= 2) ? (ydp - ydm) : 1.0;
    double domain_z = (rank >= 3) ? (zdp - zdm) : 1.0;
    double box_size = domain_x;
    double box_mass = domain_x * domain_y * domain_z * density_initial_;

    float v_rms = mach_number_;

    edot_ = 0.81/box_size*box_mass*v_rms*v_rms*v_rms;
    edot_  *= 0.8;
  }

  double norm = 0.0;

  if (edot_ != 0.0) {

    double vad = g[id_vad];
    double aad = g[id_aad];
    double zon = g[id_zones];

    const bool small_g0 = std::abs(vad) < 1e-30;

    norm = small_g0 ? 0.0001 : 1.25*dt*edot_*zon/vad;

  }

  // ASSUMES CONSTANT TIME STEP

  // double dt0 = dt;
  // norm = (dt/dt0)*norm;


  if (block->index().is_root()) {

    Monitor * monitor = cello::monitor();
    
    monitor->print ("Method","sum v*a*d    " "%.17g", g[id_vad]);
    monitor->print ("Method","sum a*a*d    " "%.17g", g[id_aad]);
    monitor->print ("Method","sum v*v*d    " "%.17g", g[id_vvd]);
    monitor->print ("Method","sum v*v      " "%.17g", g[id_vv]);
    monitor->print ("Method","sum b*b      " "%.17g", g[id_bb]);
    monitor->print ("Method","sum b*b/d    " "%.17g", g[id_bbod]);
    monitor->print ("Method","sum pr       " "%.17g", g[id_pr]);
    monitor->print ("Method","sum pr/d     " "%.17g", g[id_prod]);

    monitor->print ("Method","sum d*ax     " "%.17g", g[id_dax]);
    monitor->print ("Method","sum d*ay     " "%.17g", g[id_day]);
    monitor->print ("Method","sum d*az     " "%.17g", g[id_daz]);

    monitor->print ("Method","sum d*vx     " "%.17g", g[id_dvx]);
    monitor->print ("Method","sum d*vy     " "%.17g", g[id_dvy]);
    monitor->print ("Method","sum d*vz     " "%.17g", g[id_dvz]);

    monitor->print ("Method","sum bx       " "%.17g", g[id_bx]);
    monitor->print ("Method","sum by       " "%.17g", g[id_by]);
    monitor->print ("Method","sum bz       " "%.17g", g[id_bz]);

    monitor->print ("Method","sum d        " "%.17g", g[id_d]);
    monitor->print ("Method","sum d*d      " "%.17g", g[id_dd]);
    monitor->print ("Method","sum ln(d)    " "%.17g", g[id_lnd]);
    monitor->print ("Method","sum d*ln(d)  " "%.17g", g[id_dlnd]);
    monitor->print ("Method","min d        " "%.17g", g[id_mind]);
    monitor->print ("Method","max d        " "%.17g", g[id_maxd]);

    monitor->print ("Method","sum zones    " "%.17g", g[id_zones]);

    monitor->print ("Method","norm         " "%.17g", norm);
    monitor->print ("Method","gamma        " "%.17g", gamma_);

    monitor->print ("Method","kinetic energy             " "%.17g",
		    0.50*g[id_vvd]/g[id_zones]);
    monitor->print ("Method","turbulent magnetic energy  " "%.17g",
		    0.50*(g[id_bb]/g[id_zones]-bnotx*bnotx));
    if (gamma_ != 1.0) {
    monitor->print ("Method","internal energy            " "%.17g",
		    g[id_pr]/g[id_zones]/(gamma-1.0);
		    }
    monitor->print ("Method","potential energy           " "%.17g",
		    g[id_dlnd]/g[id_zones]);
    monitor->print ("Method","zones                      " "%.17g",
		    g[id_zones]);
    monitor->print ("Method","bnotx                      " "%.17g",
		    bnotx);
    monitor->print ("Method","<d>                        " "%.17g",
		    g[id_d]/g[id_zones]);
    monitor->print ("Method","<ln(d)>                    " "%.17g",
		    g[id_lnd]/g[id_zones]);
    monitor->print ("Method","volume-weighed rms Mach_s  " "%.17g",
		    sqrt(g[id_vv]/g[id_zones]));
    monitor->print ("Method","volume-weighed rms Mach_a  " "%.17g",
		    sqrt(g[id_vv] /
		         g[id_bbod]));
    monitor->print ("Method","mass-weighted rms Mach_s   " "%.17g",
		    sqrt(g[id_vvd]/g[id_zones]));
    monitor->print ("Method","density variance           " "%.17g",
		    sqrt(g[id_dd]/g[id_zones]));
    monitor->print ("Method","<|div(b)|>                 " "%.17g",
		    g[id_divb]/8.0/g[id_zones]);
    monitor->print ("Method","density contrast           " "%.17g",
		    g[id_maxd] /
		    g[id_mind]);                  
  }

  if (block->is_leaf()) {
    compute_resume_(block,msg);
  }

  enzo_block->compute_done();

}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceMhdIG::compute_resume_ 
(Block * block, CkReductionMsg * msg) throw()
{
  double * g = (double *)msg->getData();
  
  TRACE_TURBULENCE;  

  // Compute normalization

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = block->data()->field();

  int mx,my,mz; // total block size
  int nx,ny,nz; // active block size
  int gx,gy,gz; // number of ghost layers on each side of the block
  field.dimensions (0,&mx,&my,&mz);
  field.size         (&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);

  double dt = block->dt();

  double vad = g[id_vad];
  double aad = g[id_aad];
  double zon = g[id_zones];

  const bool small_g0 = std::abs(vad) < 1e-30;

  double norm = small_g0 ? 0.0001 : 1.25*dt*edot_*zon/vad;
  
  // ASSUMES CONSTANT TIME STEP

  double dt0 = dt;
  norm = (dt/dt0)*norm;

  const int rank = (my == 1) ? 1 : ((mz == 1) ? 2 : 3);

  enzo_float * v3[3] = {
    (enzo_float*) field.values ("velox"),
    (enzo_float*) field.values ("veloy"),
    (enzo_float*) field.values ("veloz") };
  enzo_float * v3_rx[3] = {
    (enzo_float*) field.values ("velox_rx"),
    (enzo_float*) field.values ("veloy_rx"),
    (enzo_float*) field.values ("veloz_rx") };
  enzo_float * v3_ry[3] = {
    (enzo_float*) field.values ("velox_ry"),
    (enzo_float*) field.values ("veloy_ry"),
    (enzo_float*) field.values ("veloz_ry") };
  enzo_float * v3_rz[3] = {
    (enzo_float*) field.values ("velox_rz"),
    (enzo_float*) field.values ("veloy_rz"),
    (enzo_float*) field.values ("veloz_rz") };
  enzo_float * a3[3] = {
    (enzo_float*) field.values ("drivx"),
    (enzo_float*) field.values ("drivy"),
    (enzo_float*) field.values ("drivz") };
  enzo_float * a3_rx[3] = {
    (enzo_float*) field.values ("drivx_rx"),
    (enzo_float*) field.values ("drivy_rx"),
    (enzo_float*) field.values ("drivz_rx") };
  enzo_float * a3_ry[3] = {
    (enzo_float*) field.values ("drivx_ry"),
    (enzo_float*) field.values ("drivy_ry"),
    (enzo_float*) field.values ("drivz_ry") };
  enzo_float * a3_rz[3] = {
    (enzo_float*) field.values ("drivx_rz"),
    (enzo_float*) field.values ("drivy_rz"),
    (enzo_float*) field.values ("drivz_rz") };

  // compute components of injected bulk momentum <d*a>

  const enzo_float bm[3] = 
    { enzo_float(g[id_dax]/g[id_zones]),
      enzo_float(g[id_day]/g[id_zones]),
      enzo_float(g[id_daz]/g[id_zones]) };

  // apply forcing
  // only active zones are updated (assuming mean density of 1)

  int mx = (rank >= 1) ? nx + 2*gx : nx;
  int my = (rank >= 2) ? ny + 2*gy : ny;
  int mz = (rank >= 3) ? nz + 2*gz : nz;
  for (int i=0; i<mx*my*mz; i++) {
    for (int id=0; id<rank; id++) {
      v3[id][i]    += (a3[id][i]-bm[id])*norm;
      v3_rx[id][i] += (a3_rx[id][i]-bm[id])*norm;
      v3_ry[id][i] += (a3_ry[id][i]-bm[id])*norm;
      v3_rz[id][i] += (a3_rz[id][i]-bm[id])*norm;
    }
  }

  /*
  for (int iz=gz; iz<gz+nz; iz++) {
    for (int iy=gy; iy<gy+ny; iy++) {
      for (int ix=gx; ix<gx+nx; ix++) {
	int i = ix + mx*(iy + my*iz);
	for (int id=0; id<rank; id++) {
	  v3[id][i]    += (a3[id][i]-bm[id])*norm;
	  v3_rx[id][i] += (a3_rx[id][i]-bm[id])*norm;
	  v3_ry[id][i] += (a3_ry[id][i]-bm[id])*norm;
	  v3_rz[id][i] += (a3_rz[id][i]-bm[id])*norm;
	}
      }
    }
  }
  */

  TRACE_TURBULENCE;

  delete msg;
  block->compute_done();
}
