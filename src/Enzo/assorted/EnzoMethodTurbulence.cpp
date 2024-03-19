// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulence.cpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:31:04 UTC 2014
/// @brief    Implements the EnzoMethodTurbulence class

#include "Enzo/assorted/assorted.hpp"

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"

// #define DEBUG_TURBULENCE

#ifdef DEBUG_TURBULENCE
#   define TRACE_TURBULENCE CkPrintf ("%s:%d TRACE DEBUG_TURBULENCE\n",__FILE__,__LINE__);
#else
#   define TRACE_TURBULENCE /*   */
#endif

//----------------------------------------------------------------------

EnzoMethodTurbulence::EnzoMethodTurbulence
(double edot,
 double density_initial,
 double temperature_initial,
 double mach_number,
 bool comoving_coordinates)
  : Method(),
    density_initial_(density_initial),
    temperature_initial_(temperature_initial),
    edot_(edot),
    mach_number_(mach_number),
    comoving_coordinates_(comoving_coordinates)
{
  TRACE_TURBULENCE;

  const int rank = cello::rank();

  cello::define_field("density");
  cello::define_field("temperature");
  cello::define_field("total_energy");
  if (rank >= 1) {
    cello::define_field("velocity_x");
    cello::define_field("acceleration_x");
  }
  if (rank >= 2) {
    cello::define_field("velocity_y");
    cello::define_field("acceleration_y");
  }
  if (rank >= 3) {
    cello::define_field("velocity_z");
    cello::define_field("acceleration_z");
  }

  // Initialize default Refresh object

  cello::simulation()->refresh_set_name(ir_post_,name());

  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  // TURBULENCE parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | edot_;
  p | density_initial_;
  p | temperature_initial_;
  p | mach_number_;
  p | comoving_coordinates_;

}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::compute ( Block * block) throw()
{
  TRACE_TURBULENCE;

  EnzoBlock * enzo_block = enzo::block(block);

  Field field = block->data()->field();

  const EnzoConfig * enzo_config = enzo::config();

  EnzoComputeTemperature compute_temperature(enzo::fluid_props(),
                                             comoving_coordinates_);

  compute_temperature.compute(enzo_block);

  enzo_float *  density = (enzo_float *) field.values("density");
  enzo_float *  velocity[3] = {
    (enzo_float *) field.values("velocity_x"),
    (enzo_float *) field.values("velocity_y"),
    (enzo_float *) field.values("velocity_z") };
  enzo_float * driving[3] = {
    (enzo_float *) field.values("driving_x"),
    (enzo_float *) field.values("driving_y"),
    (enzo_float *) field.values("driving_z") };
  enzo_float * temperature = (enzo_float *) field.values("temperature");

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);
  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;

  const int n = max_turbulence_array;
  double g[n];

  for (int i=0; i<max_turbulence_array-2; i++) g[i] = 0.0;

  g[index_turbulence_mind] = std::numeric_limits<double>::max();
  g[index_turbulence_maxd] = - std::numeric_limits<double>::max();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  EnzoUnits* enzo_units = enzo::units();
  const enzo_float kelvin_per_energy_u = enzo_units->kelvin_per_energy_units();

  if (block->is_leaf()) {

    for (int iz=0; iz<nz; iz++) {
      for (int iy=0; iy<ny; iy++) {
	for (int ix=0; ix<nx; ix++) {

	  int i = (ix+gx) + ndx*((iy+gy) + ndy*(iz+gz));

	  enzo_float d  = density[i];
	  for (int id=0; id<rank; id++) {
	    enzo_float v  = velocity[id][i];
	    enzo_float v2 = v*v;
	    enzo_float a  = driving[id][i];
	    enzo_float ti = kelvin_per_energy_u  / temperature[i];

	    g[index_turbulence_vad] +=   v*a*d;
	    g[index_turbulence_aad] +=   a*a*d;
	    g[index_turbulence_vvdot] += v2*d*ti;
	    g[index_turbulence_vvot] +=  v2*ti;
	    g[index_turbulence_vvd] +=   v2*d;
	    g[index_turbulence_vv] +=    v2;
	  }
	  g[index_turbulence_dd]  +=   d*d;
	  g[index_turbulence_d]   +=   d;
	  g[index_turbulence_dax] +=  d*driving[0][i];
	  g[index_turbulence_day] +=  (rank >= 2) ? d*driving[1][i] : 0.0;
	  g[index_turbulence_daz] +=  (rank >= 3) ? d*driving[2][i] : 0.0;
	  g[index_turbulence_dvx] +=  d*velocity[0][i];
	  g[index_turbulence_dvy] +=  (rank >= 2) ? d*velocity[1][i] : 0.0;
	  g[index_turbulence_dvz] +=  (rank >= 3) ? d*velocity[2][i] : 0.0;
	  g[index_turbulence_dlnd] += d*log(d);
	  g[index_turbulence_zones] += 1;
	  g[index_turbulence_mind] =
	    std::min(g[index_turbulence_mind], (double) d);
	  g[index_turbulence_maxd] =
	    std::max(g[index_turbulence_maxd], (double) d);
	}
      }
    }
  }

  CkCallback callback (CkIndex_EnzoBlock::r_method_turbulence_end(NULL),
		       enzo_block->proxy_array());
  enzo_block->contribute(n*sizeof(double),g,r_method_turbulence_type,callback);
}

//----------------------------------------------------------------------

CkReduction::reducerType r_method_turbulence_type;

void register_method_turbulence(void)
{ r_method_turbulence_type = CkReduction::addReducer(r_method_turbulence); }

CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs)
{
  double accum[max_turbulence_array];
  for (int i=0; i<max_turbulence_array; i++) {
    accum[i] = 0.0;
  }
  accum[index_turbulence_mind] = std::numeric_limits<double>::max();
  accum[index_turbulence_maxd] = - std::numeric_limits<double>::max();

  for (int i=0; i<n; i++) {
    double * values = (double *) msgs[i]->getData();
    for (int ig=0; ig<max_turbulence_array-2; ig++) {
      accum [ig] += values[ig];
    }
    accum [index_turbulence_mind] =
      std::min(accum[index_turbulence_mind],values[index_turbulence_mind]);
    accum [index_turbulence_maxd] =
      std::max(accum[index_turbulence_maxd],values[index_turbulence_maxd]);
  }
  return CkReductionMsg::buildNew(max_turbulence_array*sizeof(double),accum);
}

//----------------------------------------------------------------------

void EnzoBlock::r_method_turbulence_end(CkReductionMsg * msg)
{
  TRACE_TURBULENCE;
  performance_start_(perf_compute,__FILE__,__LINE__);
  method()->compute_resume (this,msg);
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::compute_resume
(Block * block,
 CkReductionMsg * msg) throw()
{
  TRACE_TURBULENCE;

  double * g = (double *)msg->getData();

  Data * data = block->data();
  Field field = data->field();


  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int n = nx*ny*nz;

  double dt = block->state()->dt();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  double xdm,ydm,zdm;
  data->lower(&xdm,&ydm,&zdm);
  double xdp,ydp,zdp;
  data->upper(&xdp,&ydp,&zdp);

  // compute edot (TurbulenceSimulationInitialize.C)

  // If RandomForcingEdot (i.e. the energy injection rate) is not set
  // in the parameter file, get it from [MacLow1999] formula.  Note:
  // the formula is calibrated for generic forcing fields; coefficient
  // 0.81 can potentially be inappropriate for a purely solenoidal
  // forcing; also our Gamma is not quite 1.0.

  if (edot_ < 0.0) {
    // Only compute if needed at the beginning--could/should be in
    // EnzoInitialTurbulence
    double domain_x =               (xdp - xdm);
    double domain_y = (rank >= 2) ? (ydp - ydm) : 1.0;
    double domain_z = (rank >= 3) ? (zdp - zdm) : 1.0;
    double box_size = domain_x;
    double box_mass = domain_x * domain_y * domain_z * density_initial_;


    EnzoUnits* enzo_units = enzo::units();
    float v_rms = mach_number_ * sqrt(temperature_initial_ /
                                      enzo_units->kelvin_per_energy_units());

    edot_ = 0.81/box_size*box_mass*v_rms*v_rms*v_rms;

    // Approximate correction to the MacLow's factor (see eqs (7) - (8))
    // for **this PPM implementation**. Seems to be OK for 64^3, 128^3
    // and 256^3 Mach=3,6,10 simulations of **solenoidally** driven
    // turbulence. */

    //
    // (7) $\dot{E}_{\textsf{\scriptsize{kin}}} \simeq - \eta_{\nu} m
    //      \tilde{k} v^{3}_{\textsf{\scriptsize{rms}}}$
    //
    //
    // (8) $\dot{E}_{\textsf{\scriptsize{kin}}} = - \eta_{e} m^{-1/2}
    //      \tilde{k} E^{3/2}_{\textsf{\scriptsize{kin}}}$
    //

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
    //  //  small push at the start, when gv0==0 due to zero initial velocities
    //   if (gv0 < 1e-30 && gv0 > -1e-30 && MetaData->TopGridRank == 2) *norm = 0.0001;
    //     else    *norm = 1.25*dt*RandomForcingEdot*numberOfGridZones/gv0;


    double vad = g[index_turbulence_vad];

    const bool small_g0 = std::abs(vad) < 1e-30;

    norm = small_g0 ? 0.0001 : 1.25*dt*edot_*n/vad;

      // OLD COMPUTATION:
      //
      //      norm = ( sqrt(g[0]*g[0] + 2.0*n*g[1]*dt*edot_) - g[0] ) / g[1];
  }

  // ASSUMES CONSTANT TIME STEP

  // double dt0 = dt;
  // norm = (dt/dt0)*norm;


  if (block->index().is_root()) {

    Monitor * monitor = cello::monitor();

    monitor->print ("Method","sum v*a*d    " "%.17g", g[index_turbulence_vad]);
    monitor->print ("Method","sum a*a*d    " "%.17g",g[index_turbulence_aad]);
    monitor->print ("Method","sum v*v*d/t  " "%.17g",g[index_turbulence_vvdot]);
    monitor->print ("Method","sum v*v/t    " "%.17g",g[index_turbulence_vvot]);
    monitor->print ("Method","sum v*v*d    " "%.17g",g[index_turbulence_vvd]);
    monitor->print ("Method","sum v*v      " "%.17g",g[index_turbulence_vv]);
    monitor->print ("Method","sum d*d      " "%.17g",g[index_turbulence_dd]);

    monitor->print ("Method","sum d*ax     " "%.17g",g[index_turbulence_dax]);
    monitor->print ("Method","sum d*ay     " "%.17g",g[index_turbulence_day]);
    monitor->print ("Method","sum d*az     " "%.17g",g[index_turbulence_daz]);

    monitor->print ("Method","sum d*vx     " "%.17g",g[index_turbulence_dvx]);
    monitor->print ("Method","sum d*vy     " "%.17g",g[index_turbulence_dvy]);
    monitor->print ("Method","sum d*vz     " "%.17g",g[index_turbulence_dvz]);

    monitor->print ("Method","sum d*ln(d)  " "%.17g",g[index_turbulence_dlnd]);
    monitor->print ("Method","sum zones    " "%.17g",g[index_turbulence_zones]);

    monitor->print ("Method","min d        " "%.17g",g[index_turbulence_mind]);
    monitor->print ("Method","max d        " "%.17g",g[index_turbulence_maxd]);
    monitor->print ("Method","sum d        " "%.17g",g[index_turbulence_d]);
    monitor->print ("Method","norm         " "%.17g",norm);

    monitor->print ("Method","kinetic energy          " "%.17g",
		    0.50*g[index_turbulence_vvd]/n);
    monitor->print ("Method","mass weighted rms Mach  " "%.17g",
		    sqrt(g[index_turbulence_vvdot]/n));
    monitor->print ("Method","volume weighed rms Mach " "%.17g",
		    sqrt(g[index_turbulence_vvot]/n));
    monitor->print ("Method","rms Velocity            " "%.17g",
		    sqrt(g[index_turbulence_vv]/n));
    monitor->print ("Method","Density variance        " "%.17g",
		    sqrt(g[index_turbulence_dd]/n));
    monitor->print ("Method","min/max Density         " "%.17g",
		    g[index_turbulence_mind] /
		    g[index_turbulence_maxd]);
  }

  if (block->is_leaf()) {
    compute_resume_(block,msg);
  }

  delete msg;
  block->compute_done();

}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::compute_resume_
(Block * block, CkReductionMsg * msg) throw()
{

  TRACE_TURBULENCE;
  // Compute normalization

  Field field = block->data()->field();

  int mx,my,mz;
  int nx,ny,nz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.size         (&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);

  int n = nx*ny*nz;

  double * g = (double *)msg->getData();

  double dt = block->state()->dt();

  double norm = (edot_ != 0.0) ?
    ( sqrt(g[0]*g[0] + 2.0*n*g[1]*dt*edot_) - g[0] ) / g[1] : 0.0;

  // ASSUMES CONSTANT TIME STEP

  double dt0 = dt;
  norm = (dt/dt0)*norm;

  const int rank = cello::rank();

  enzo_float * te = (enzo_float*) field.values ("total_energy");
  enzo_float * v3[3] = {
    (enzo_float*) field.values ("velocity_x"),
    (enzo_float*) field.values ("velocity_y"),
    (enzo_float*) field.values ("velocity_z") };
  enzo_float * a3[3] = {
    (enzo_float*) field.values ("driving_x"),
    (enzo_float*) field.values ("driving_y"),
    (enzo_float*) field.values ("driving_z") };

  // compute bulk momentum
  const enzo_float bm[3] =
    { enzo_float(g[index_turbulence_dax]/n),
      enzo_float(g[index_turbulence_day]/n),
      enzo_float(g[index_turbulence_daz]/n)};

  //  for (int dim = 0; dim <  MetaData->TopGridRank; dim++)
  //	bulkMomentum[dim] = GlobVal[7+dim]/numberOfGridZones;

  for (int iz=gz; iz<gz+nz; iz++) {
    for (int iy=gy; iy<gy+ny; iy++) {
      for (int ix=gx; ix<gx+nx; ix++) {
	int i = ix + mx*(iy + my*iz);
	for (int id=0; id<rank; id++) {
	  //	  	  te[i] += v3[id][i]*a3[id][i]*norm + 0.5*a3[id][i]*norm*a3[id][i]*norm;
	  //	  	  v3[id][i] += a3[id][i]*norm;
	  te[i] += (v3[id][i]*(a3[id][i]-bm[id]))*norm;
	  v3[id][i] += (a3[id][i]-bm[id])*norm;

	}
      }
    }
  }

  TRACE_TURBULENCE;

}
