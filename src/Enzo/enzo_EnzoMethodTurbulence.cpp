// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulence.cpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:31:04 UTC 2014
/// @brief    Implements the EnzoMethodTurbulence class

#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

// #define DEBUG_TURBULENCE

#ifdef DEBUG_TURBULENCE
#   define TRACE_TURBULENCE CkPrintf ("%s:%d TRACE DEBUG_TURBULENCE\n",__FILE__,__LINE__);
#else
#   define TRACE_TURBULENCE /*   */
#endif

//----------------------------------------------------------------------

EnzoMethodTurbulence::EnzoMethodTurbulence 
(const FieldDescr * field_descr,
 double edot,
 double density_initial,
 double temperature_initial,
 double mach_number,
 int comoving_coordinates)
  : Method(),
    density_initial_(density_initial),
    temperature_initial_(temperature_initial),
    edot_(edot),
    mach_number_(mach_number),
    comoving_coordinates_(comoving_coordinates)
{
  TRACE_TURBULENCE;  
  
  // Initialize default Refresh object

  const int ir = add_refresh(4,0,neighbor_leaf,sync_barrier);

  refresh(ir)->add_all_fields(field_descr->field_count());

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

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);
  
  Field field = block->data()->field();

  const EnzoConfig * enzo_config = static_cast<const EnzoConfig*>
    (enzo_block->simulation()->config());

  EnzoComputeTemperature compute_temperature
    (enzo_config->ppm_density_floor,
     enzo_config->ppm_temperature_floor,
     enzo_config->ppm_mol_weight,
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

  const int n = sizeof(enzo_block->method_turbulence_data)/sizeof(double);
  double * g = enzo_block->method_turbulence_data;
  
  ASSERT2 ("EnzoMethodTurbulence::compute()",
	   "Size of EnzoBlock::method_turbulence_data array %d "
	   "must be at least %d",
	   n, max_turbulence_array, 
	   (n >= max_turbulence_array));

  for (int i=0; i<max_turbulence_array-2; i++) g[i] = 0.0;

  g[index_turbulence_mind] = std::numeric_limits<double>::max();
  g[index_turbulence_maxd] = - std::numeric_limits<double>::max();

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

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
	    enzo_float ti = 1.0 / temperature[i];

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
	  g[index_turbulence_mind] =
	    std::min(g[index_turbulence_mind], (double) d);
	  g[index_turbulence_maxd] =
	    std::max(g[index_turbulence_maxd], (double) d);
	}
      }
    }
  }
  enzo_block->method_turbulence_begin();
}

//----------------------------------------------------------------------

extern CkReduction::reducerType r_method_turbulence_type;

void EnzoBlock::method_turbulence_begin()
{
  TRACE_TURBULENCE;  
  const int n = max_turbulence_array * sizeof(double);
  double * g = method_turbulence_data;
  CkCallback callback (CkIndex_EnzoBlock::p_method_turbulence_end(NULL),
		       thisProxy);
  contribute(n,g,r_method_turbulence_type,callback);
}

//----------------------------------------------------------------------

// SEE main.cpp for implementation

// CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs)

//----------------------------------------------------------------------

void EnzoBlock::p_method_turbulence_end(CkReductionMsg * msg)
{
  TRACE_TURBULENCE;  
  performance_start_(perf_compute,__FILE__,__LINE__);
  method()->compute_resume (this,msg);
  delete msg;
  performance_stop_(perf_compute,__FILE__,__LINE__);
}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::compute_resume 
(Block * block,
 CkReductionMsg * msg) throw()
{
  TRACE_TURBULENCE;  

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  double * g = enzo_block->method_turbulence_data;
  for (int i=0; i<max_turbulence_array; i++) {
    g[i] = ((double *)msg->getData())[i];
  }

  Data * data = block->data();
  Field field = data->field();
  

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int n = nx*ny*nz;

  double dt = block->dt();

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

    float v_rms = mach_number_ * sqrt(temperature_initial_);

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

    Monitor * monitor = block->simulation()->monitor();
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

  if (!block->is_leaf()) {
    enzo_block->compute_done();
    return;
  }

  const int p = field.precision (0);

  if      (p == precision_single)    
    compute_resume_<float> (block,msg);
  else if (p == precision_double)    
    compute_resume_<double> (block,msg);
  else if (p == precision_quadruple) 
    compute_resume_<long double> (block,msg);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodTurbulence::compute_resume_ 
(Block * block,
 CkReductionMsg * msg) throw()
{

  TRACE_TURBULENCE;  
  // Compute normalization

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = block->data()->field();

  int ndx,ndy,ndz;
  field.size(&ndx,&ndy,&ndz);
  int nd = ndx*ndy*ndz;

  double * g = enzo_block->method_turbulence_data;

  double dt = block->dt();

  double norm = (edot_ != 0.0) ?
    ( sqrt(g[0]*g[0] + 2.0*nd*g[1]*dt*edot_) - g[0] ) / g[1] : 0.0;

  // ASSUMES CONSTANT TIME STEP

  double dt0 = dt;
  norm = (dt/dt0)*norm;

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);
  const int rank = ((mz == 1) ? ((my == 1) ? 1 : 2) : 3);

  T * te = (T*) field.values ("total_energy");
  T * v3[3] = { (T*) field.values ("velocity_x"),
		(T*) field.values ("velocity_y"),
		(T*) field.values ("velocity_z") };
  T * a3[3] = { (T*) field.values ("driving_x"),
		(T*) field.values ("driving_y"),
		(T*) field.values ("driving_z") };

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);
  int gx,gy,gz;
  field.ghost_depth(0,&gx,&gy,&gz);
  int nbx = nx + 2*gx;
  int nby = ny + 2*gy;

  // compute bulk momentum
  const T bm[3] = 
    { T(g[index_turbulence_dax]/nd),
      T(g[index_turbulence_day]/nd),
      T(g[index_turbulence_daz]/nd)};

  //  for (int dim = 0; dim <  MetaData->TopGridRank; dim++)
  //	bulkMomentum[dim] = GlobVal[7+dim]/numberOfGridZones;

  for (int iz=gz; iz<gz+nz; iz++) {
    for (int iy=gy; iy<gy+ny; iy++) {
      for (int ix=gx; ix<gx+nx; ix++) {
	int i = ix + nbx*(iy + nby*iz);
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
  enzo_block->compute_done();
  
}
