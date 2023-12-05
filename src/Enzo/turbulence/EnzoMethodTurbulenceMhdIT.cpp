// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulenceMhdIT.cpp
/// @author   Alexei Kritsuk (akritsuk@ucsd.edu)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:31:04 UTC 2014
/// @date     Fri Aug 24 00:31:04 UTC 2018
/// @brief    Implements the EnzoMethodTurbulenceMhdIT class

#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

//----------------------------------------------------------------------

enum {
  it_vad,
  it_aad,
  it_vvd,
  it_vv,
  it_dvx,
  it_dvy,
  it_dvz,
  it_dax,
  it_day,
  it_daz,
  it_bx,
  it_by,
  it_bz,
  it_bb,
  it_bbod,
  it_divb,
  it_d,
  it_dd,
  it_lnd,
  it_dlnd,
  it_zones,
  it_mind,
  it_maxd,
  max_turbulence_mhd_it_array
};

//#define DEBUG_TURBULENCE

#ifdef DEBUG_TURBULENCE
#   define TRACE_TURBULENCE CkPrintf ("%s:%d TRACE DEBUG_TURBULENCE\n",__FILE__,__LINE__);
#else
#   define TRACE_TURBULENCE /*   */
#endif

//----------------------------------------------------------------------

EnzoMethodTurbulenceMhdIT::EnzoMethodTurbulenceMhdIT 
(double edot,
 double density_initial,
 double bfieldx_initial,
 double mach_number,
 bool comoving_coordinates)
  : Method(),
    density_initial_(density_initial),
    bfieldx_initial_(bfieldx_initial),
    edot_(edot),
    mach_number_(mach_number),
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

void EnzoMethodTurbulenceMhdIT::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | edot_;
  p | density_initial_;
  p | bfieldx_initial_;
  p | mach_number_;
  p | comoving_coordinates_;

}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceMhdIT::compute ( Block * block) throw()
{
  TRACE_TURBULENCE;  

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  
  Field field = block->data()->field();

  enzo_float *  density = (enzo_float *) field.values("density");
  enzo_float *  dens_rx = (enzo_float *) field.values("dens_rx");
  enzo_float *  dens_ry = (enzo_float *) field.values("dens_ry");
  enzo_float *  dens_rz = (enzo_float *) field.values("dens_rz");

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
  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;

  const int n = max_turbulence_mhd_it_array;
  double * g = new double[n];

  for (int i=0; i<max_turbulence_mhd_it_array-2; i++) g[i] = 0.0;

  g[it_mind] =   std::numeric_limits<double>::max();
  g[it_maxd] = - std::numeric_limits<double>::max();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  if (block->is_leaf()) {

    // Loops cover only active zones of this block and
    // compute averages for zone centers (mostly) to be used for forcing normalization

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {

	  int i = (ix+gx) + ndx*((iy+gy) + ndy*(iz+gz));

	  enzo_float d   = density[i];
	  for (int id=0; id<rank; id++) {
	    enzo_float v  = velocity[id][i];
	    enzo_float v2 = v*v;
	    enzo_float a  = driving[id][i];
	    enzo_float b  = bfield[id][i];
	    enzo_float b2 = b*b;

	    g[it_vad]   += v*a*d; // 0
	    g[it_aad]   += a*a*d; // 1
	    g[it_vvd]   += v2*d;  // 2
	    g[it_vv]    += v2;    // 3
	    g[it_bb]    += b2;    // 13
	    g[it_bbod]  += b2/d;  // 14
	  }
	  g[it_d]    +=  d;       // 16
	  g[it_dd]   +=  d*d;     // 17
	  g[it_lnd]  +=  log(d);  // 18
	  g[it_dlnd] +=  d*log(d);// 19

          g[it_zones] += 1;       // 20

	  g[it_dvx]  +=  d*velocity[0][i];                     // 4
	  g[it_dvy]  +=  (rank >= 2) ? d*velocity[1][i] : 0.0; // 5
	  g[it_dvz]  +=  (rank >= 3) ? d*velocity[2][i] : 0.0; // 6
	  /*
            g[it_dax]  +=  d*driving[0][i];                      // 7
            g[it_day]  +=  (rank >= 2) ? d*driving[1][i] : 0.0;  // 8
            g[it_daz]  +=  (rank >= 3) ? d*driving[2][i] : 0.0;  // 9
	  */
	  g[it_dax]  +=  ( d*driving[0][i] + dens_rx[i]*driv_rx[0][i] + 
                                                dens_ry[i]*driv_ry[0][i] + 
                                                dens_rz[i]*driv_rz[0][i] )/4.0;                      // 7
	  g[it_day]  +=  (rank >= 2) ? ( d*driving[1][i] + dens_rx[i]*driv_rx[1][i] +
                                                              dens_ry[i]*driv_ry[1][i] +
                                                              dens_rz[i]*driv_rz[1][i] )/4.0 : 0.0;  // 8
	  g[it_daz]  +=  (rank >= 3) ? ( d*driving[2][i] + dens_rx[i]*driv_rx[2][i] +
                                                              dens_ry[i]*driv_ry[2][i] +
                                                              dens_rz[i]*driv_rz[2][i] )/4.0 : 0.0;  // 9

	  g[it_bx]   +=  bfield[0][i];                         // 10
	  g[it_by]   +=  (rank >= 2) ? bfield[1][i] : 0.0;     // 11
	  g[it_bz]   +=  (rank >= 3) ? bfield[2][i] : 0.0;     // 12

	  // PPML-style divergence calculation in 3D.
	  // One can also use max|div(b)| for control, see Fig. 18 in Ustyugov et al. (2009, JCP 228, 7614).
	  int is = 1;
	  int js = ndx;
	  int ks = ndx*ndy;
	  g[it_divb] +=  fabs(bfield_rx[0][i+is] - bfield_rx[0][i-is] +          // 15
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

	  g[it_mind] =                                         // 21
	    std::min(g[it_mind], (double) d);
	  g[it_maxd] =                                         // 22
	    std::max(g[it_maxd], (double) d);
	}
      }
    }
  }

  CkCallback callback (CkIndex_EnzoBlock::r_method_turbulence_mhd_it_end(NULL),
                       enzo_block->proxy_array());
  enzo_block->contribute(n*sizeof(double),g,r_method_turbulence_mhd_it_type,callback);

  delete [] g;

}

//----------------------------------------------------------------------

CkReduction::reducerType r_method_turbulence_mhd_it_type;

void register_method_turbulence_mhd_it(void)
{ 
  r_method_turbulence_mhd_it_type = CkReduction::addReducer(r_method_turbulence_mhd_it); 
}

CkReductionMsg * r_method_turbulence_mhd_it(int n, CkReductionMsg ** msgs)
{
  double accum[max_turbulence_mhd_it_array];
  for (int i=0; i<max_turbulence_mhd_it_array; i++) {
    accum[i] = 0.0;
  }
  accum[it_mind] =   std::numeric_limits<double>::max();
  accum[it_maxd] = - std::numeric_limits<double>::max();

  for (int i=0; i<n; i++) {
    double * values = (double *) msgs[i]->getData();
    for (int ig=0; ig<max_turbulence_mhd_it_array-2; ig++) {
      accum [ig] += values[ig];
    }
    accum [it_mind] = 
      std::min(accum[it_mind],values[it_mind]);
    accum [it_maxd] = 
      std::max(accum[it_maxd],values[it_maxd]);
  }
  return CkReductionMsg::buildNew(max_turbulence_mhd_it_array*sizeof(double),accum);
}


//----------------------------------------------------------------------

void EnzoBlock::r_method_turbulence_mhd_it_end(CkReductionMsg * msg)
{
  TRACE_TURBULENCE;  
  performance_start_(perf_compute,__FILE__,__LINE__);
  method()->compute_resume (this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceMhdIT::compute_resume 
(Block * block,
 CkReductionMsg * msg) throw()
{
  TRACE_TURBULENCE;  

  double * g = (double *)msg->getData();

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

  // compute edot (TurbulenceSimulationInitialize.C)

  // If RandomForcingEdot (i.e. the energy injection rate) is not set
  // in the parameter file, get it from the MacLow (1999) formula.  Note:
  // the formula is calibrated for generic forcing fields; coefficient
  // 0.81 can potentially be inappropriate for a purely solenoidal
  // forcing.

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
 
    // Approximate correction to the MacLow's factor (see eqs (7) - (8))
    // for Enzo's PPM implementation. Seems to be OK for 64^3, 128^3
    // and 256^3 Mach=3,6,10 simulations of solenoidally driven
    // turbulence.
    //
    // (7) $\dot{E}_{\textsf{\scriptsize{kin}}} \simeq - \eta_{\nu} m
    //      \tilde{k} v^{3}_{\textsf{\scriptsize{rms}}}$
    //
    //
    // (8) $\dot{E}_{\textsf{\scriptsize{kin}}} = - \eta_{e} m^{-1/2}
    //      \tilde{k} E^{3/2}_{\textsf{\scriptsize{kin}}}$
 
    edot_  *= 0.8;
 
  }

  // compute norm (ComputeRandomForcingNormalization.C)

  double norm = 0.0;

  if (edot_ != 0.0) {

    // Original code in ComputeRandomForcingNormalization.C:
    //
    //   float gv0 = GlobVal[0];
    //   if (gv0 < 1e-30 && gv0 > -1e-30 && MetaData->TopGridRank == 3) {ERROR_MESSAGE} 
    //      else    *norm = 1.25*dt*RandomForcingEdot*numberOfGridZones/gv0;
    //      small push at the start, when gv0==0 due to zero initial velocities
    //   if (gv0 < 1e-30 && gv0 > -1e-30 && MetaData->TopGridRank == 2) *norm = 0.0001;
    //     else    *norm = 1.25*dt*RandomForcingEdot*numberOfGridZones/gv0;

    
    double vad = g[it_vad];
    double aad = g[it_aad];
    double zon = g[it_zones];

    const bool small_g0 = std::abs(vad) < 1e-30;

    norm = small_g0 ? 0.0001 : 1.25*dt*edot_*zon/vad;

    //    norm = (edot_ != 0.0) ? (sqrt(vad*vad + 2.0*zon*aad*dt*edot_) - vad)/aad : 0.0;
   
      
    // OLD COMPUTATION:
    //
    //      norm = ( sqrt(g[0]*g[0] + 2.0*n*g[1]*dt*edot_) - g[0] ) / g[1];
  }

  // ASSUMES CONSTANT TIME STEP

  // double dt0 = dt;
  // norm = (dt/dt0)*norm;

  monitor_output_(block,g,norm,bnotx);

  if (block->is_leaf()) {
    compute_resume_(block,msg);
  }

  delete msg;
  block->compute_done();

}

//======================================================================

void EnzoMethodTurbulenceMhdIT::monitor_output_
(Block * block, double * g, double norm, double bnotx)
{
  if (block->index().is_root()) {

    Monitor * monitor = cello::monitor();
    
    monitor->print ("Method","sum v*a*d    " "%.17g", g[it_vad]);
    monitor->print ("Method","sum a*a*d    " "%.17g", g[it_aad]);
    monitor->print ("Method","sum v*v*d    " "%.17g", g[it_vvd]);
    monitor->print ("Method","sum v*v      " "%.17g", g[it_vv]);
    monitor->print ("Method","sum b*b      " "%.17g", g[it_bb]);
    monitor->print ("Method","sum b*b/d    " "%.17g", g[it_bbod]);

    monitor->print ("Method","sum d*ax     " "%.17g", g[it_dax]);
    monitor->print ("Method","sum d*ay     " "%.17g", g[it_day]);
    monitor->print ("Method","sum d*az     " "%.17g", g[it_daz]);

    monitor->print ("Method","sum d*vx     " "%.17g", g[it_dvx]);
    monitor->print ("Method","sum d*vy     " "%.17g", g[it_dvy]);
    monitor->print ("Method","sum d*vz     " "%.17g", g[it_dvz]);

    monitor->print ("Method","sum bx       " "%.17g", g[it_bx]);
    monitor->print ("Method","sum by       " "%.17g", g[it_by]);
    monitor->print ("Method","sum bz       " "%.17g", g[it_bz]);

    monitor->print ("Method","sum d        " "%.17g", g[it_d]);
    monitor->print ("Method","sum d*d      " "%.17g", g[it_dd]);
    monitor->print ("Method","sum ln(d)    " "%.17g", g[it_lnd]);
    monitor->print ("Method","sum d*ln(d)  " "%.17g", g[it_dlnd]);
    monitor->print ("Method","min d        " "%.17g", g[it_mind]);
    monitor->print ("Method","max d        " "%.17g", g[it_maxd]);

    monitor->print ("Method","sum zones    " "%.17g", g[it_zones]);

    monitor->print ("Method","norm         " "%.17g", norm);
    monitor->print ("Method","edot         " "%.17g", edot_);

    monitor->print ("Method","kinetic energy             " "%.17g",
		    0.50*g[it_vvd]/g[it_zones]);
    monitor->print ("Method","turbulent magnetic energy  " "%.17g",
		    0.50*(g[it_bb]/g[it_zones]-bnotx*bnotx));
    monitor->print ("Method","potential energy           " "%.17g",
		    g[it_dlnd]/g[it_zones]);
    monitor->print ("Method","zones                      " "%.17g",
		    g[it_zones]);
    monitor->print ("Method","bnotx                      " "%.17g",
		    bnotx);
    monitor->print ("Method","<d>                        " "%.17g",
		    g[it_d]/g[it_zones]);
    monitor->print ("Method","<ln(d)>                    " "%.17g",
		    g[it_lnd]/g[it_zones]);
    monitor->print ("Method","volume-weighed rms Mach_s  " "%.17g",
		    sqrt(g[it_vv]/g[it_zones]));
    monitor->print ("Method","volume-weighed rms Mach_a  " "%.17g",
		    sqrt(g[it_vv] /
		         g[it_bbod]));
    monitor->print ("Method","mass-weighted rms Mach_s   " "%.17g",
		    sqrt(g[it_vvd]/g[it_zones]));
    monitor->print ("Method","density variance           " "%.17g",
		    sqrt(g[it_dd]/g[it_zones]));
    monitor->print ("Method","<|div(b)|>                 " "%.17g",
		    g[it_divb]/8.0/g[it_zones]);
    monitor->print ("Method","density contrast           " "%.17g",
		    g[it_maxd] /
		    g[it_mind]);                  
  }
}

//----------------------------------------------------------------------

void EnzoMethodTurbulenceMhdIT::compute_resume_ 
(Block * block, CkReductionMsg * msg) throw()
{
  
  TRACE_TURBULENCE;  

  // Compute normalization

  Field field = block->data()->field();

  int mx,my,mz; // total block size
  int nx,ny,nz; // active block size
  int gx,gy,gz; // number of ghost layers on each side of the block
  field.dimensions (0,&mx,&my,&mz);
  field.size         (&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);

  int n = nx*ny*nz;

  double * g = (double *)msg->getData();

  double dt = block->dt();

  double vad = g[it_vad];
  double aad = g[it_aad];
  double zon = g[it_zones];

  const bool small_g0 = std::abs(vad) < 1e-30;

  double norm = small_g0 ? 0.0001 : 1.25*dt*edot_*zon/vad;
  
  //  double norm = (edot_ != 0.0) ? (sqrt(vad*vad + 2.0*zon*aad*dt*edot_) - vad)/aad : 0.0;
  
  // ASSUMES CONSTANT TIME STEP

  double dt0 = dt;
  norm = (dt/dt0)*norm;

  const int rank = cello::rank();
  //  const int rank = (my == 1) ? 1 : ((mz == 1) ? 2 : 3);

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

  // compute injected bulk momentum <d*a> in x, y, and z directions

  const enzo_float bm[3] = 
    { enzo_float(g[it_dax]/g[it_zones]),
      enzo_float(g[it_day]/g[it_zones]),
      enzo_float(g[it_daz]/g[it_zones]) };

  //  if (block->index().is_root()) {
  //    Monitor * monitor = cello::monitor();
  //    monitor->print ("Method","bulk momentum    " "%.17g %.17g %.17g", bm[0], bm[1], bm[2]);
  //  }

  // apply forcing
  // only active zones are updated (assuming mean density of 1)

  int ndx = (rank >= 1) ? nx + 2*gx : nx;
  int ndy = (rank >= 1) ? ny + 2*gy : ny;
  int ndz = (rank >= 1) ? nz + 2*gz : nz;
  for (int i=0; i<ndx*ndy*ndz; i++) {
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
  
}
