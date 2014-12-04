// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulence.cpp
/// @author   Alexei Kritsuk (kritsuk@gmail.com)
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:31:04 UTC 2014
/// @brief    Implements the EnzoMethodTurbulence class

#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

//----------------------------------------------------------------------

EnzoMethodTurbulence::EnzoMethodTurbulence 
(double edot,
 double density_initial,
 double temperature_initial,
 double mach_number)
  : Method(),
    density_initial_(density_initial),
    temperature_initial_(temperature_initial),
    edot_(edot),
    mach_number_(mach_number)
{
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

}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::compute ( CommBlock * comm_block) throw()
{

  // MUST PROCEED EVEN IF NOT A LEAF NODE SINCE REDUCTION IS ON ALL
  // CommBlock's
  //  if (!comm_block->is_leaf()) return;

  // Initialize Method parameters: used for domain extents
  initialize_(comm_block);

  //  INCOMPLETE("EnzoMethodTurbulence::compute()");

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);
  
  Field field = comm_block->block()->field();

  const EnzoConfig * enzo_config = static_cast<const EnzoConfig*>
    (enzo_block->simulation()->config());

  EnzoComputeTemperature compute_temperature
    (enzo_config->ppm_density_floor,
     enzo_config->ppm_temperature_floor,
     enzo_config->ppm_mol_weight);

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
  field.ghosts(0,&gx,&gy,&gz);
  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;

  const int n = sizeof(enzo_block->method_turbulence_data)/sizeof(double);
  double * g = enzo_block->method_turbulence_data;
  
  ASSERT2 ("EnzoMethodTurbulence::compute()",
	   "Size of EnzoBlock::method_turbulence_data array %d "
	   "must be at least %d",
	   n, MAX_TURBULENCE_ARRAY, 
	   (n >= MAX_TURBULENCE_ARRAY));

  for (int i=0; i<MAX_TURBULENCE_ARRAY-2; i++) g[i] = 0.0;

  g[INDEX_TURBULENCE_minD] = std::numeric_limits<double>::max();
  g[INDEX_TURBULENCE_maxD] = std::numeric_limits<double>::min();

  const int rank = this->rank();

  if (comm_block->is_leaf()) {

    for (int iz=gz; iz<gz+nz; iz++) {
      for (int iy=gy; iy<gy+ny; iy++) {
	for (int ix=gx; ix<gx+nx; ix++) {

	  int i = ix + ndx*(iy + ndy*iz);

	  enzo_float d  = density[i];
	  for (int id=0; id<rank; id++) {
	    enzo_float v  = velocity[id][i];
	    enzo_float v2 = v*v;
	    enzo_float a  = driving[id][i];
	    enzo_float ti = 1.0 / temperature[i];

	    g[INDEX_TURBULENCE_VAD] +=   v*a*d;
	    g[INDEX_TURBULENCE_AAD] +=   a*a*d;
	    g[INDEX_TURBULENCE_VVDoT] += v2*d*ti;
	    g[INDEX_TURBULENCE_VVoT] +=  v2*ti;
	    g[INDEX_TURBULENCE_VVD] +=   v2*d;
	    g[INDEX_TURBULENCE_VV] +=    v2;
	  }
	  g[INDEX_TURBULENCE_DD] +=   d*d;
	  g[INDEX_TURBULENCE_DAx] +=  d*driving[0][i];
	  g[INDEX_TURBULENCE_DAy] +=  (rank >= 2) ? d*driving[1][i] : 0.0;
	  g[INDEX_TURBULENCE_DAz] +=  (rank >= 3) ? d*driving[2][i] : 0.0;
	  g[INDEX_TURBULENCE_DVx] +=  d*velocity[0][i];
	  g[INDEX_TURBULENCE_DVy] +=  (rank >= 2) ? d*velocity[1][i] : 0.0;
	  g[INDEX_TURBULENCE_DVz] +=  (rank >= 3) ? d*velocity[2][i] : 0.0;
	  g[INDEX_TURBULENCE_DlnD] += d*log(d);
	  g[INDEX_TURBULENCE_minD] =
	    std::min(g[INDEX_TURBULENCE_minD], (double) d);
	  g[INDEX_TURBULENCE_maxD] =
	    std::max(g[INDEX_TURBULENCE_maxD], (double) d);
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
  const int n = MAX_TURBULENCE_ARRAY * sizeof(double);
  double * g = method_turbulence_data;
  CkCallback cb (CkIndex_EnzoBlock::p_method_turbulence_end(NULL),thisProxy);
  contribute(n,g,r_method_turbulence_type,cb);
}

//----------------------------------------------------------------------

// SEE main.cpp for implementation

// CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs)

//----------------------------------------------------------------------

void EnzoBlock::p_method_turbulence_end(CkReductionMsg * msg)
{
  method()->compute_resume (this,msg);
  delete msg;
}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::compute_resume 
(CommBlock * comm_block,
 CkReductionMsg * msg) throw()
{

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);

  double * g = enzo_block->method_turbulence_data;
  for (int i=0; i<MAX_TURBULENCE_ARRAY; i++) {
    g[i] = ((double *)msg->getData())[i];
  }

  int nx,ny,nz;
  domain_size(&nx,&ny,&nz);
  int n = nx*ny*nz;

  double dt = comm_block->dt();

  const int rank = this->rank();

  double xdm,ydm,zdm;
  lower_domain(&xdm,&ydm,&zdm);
  double xdp,ydp,zdp;
  upper_domain(&xdp,&ydp,&zdp);

  // compute edot (TurbulenceSimulationInitialize.C)

  // If RandomForcingEdot (i.e. the energy injection rate) is not set
  // in the parameter file, get it from [MacLow1999] formula.  Note:
  // the formula is calibrated for generic forcing fields; coefficient
  // 0.81 can potentially be inappropriate for a purely solenoidal
  // forcing; also our Gamma is not quite 1.0.

  if (edot_ < 0.0) {
    // Only compute if needed at the beginning--could/should be in
    // EnzoInitialTurbulence
    double domain_x = (rank >= 1) ? (xdp - xdm) : 1.0;
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

    
    double vad = g[INDEX_TURBULENCE_VAD];

    const bool small_g0 = std::abs(vad < 1e-30);

    norm = small_g0 ? 0.0001 : 1.25*dt*edot_*n/vad;
      
      // OLD COMPUTATION:
      //
      //      norm = ( sqrt(g[0]*g[0] + 2.0*n*g[1]*dt*edot_) - g[0] ) / g[1];
  }

  // ASSUMES CONSTANT TIME STEP

  // double dt0 = dt;
  // norm = (dt/dt0)*norm;


  if (comm_block->index().is_root()) {

    Monitor * monitor = comm_block->simulation()->monitor();

    monitor->print ("Method","sum v*a*d   = %lg",g[INDEX_TURBULENCE_VAD]);
    monitor->print ("Method","sum a*a*d   = %lg",g[INDEX_TURBULENCE_AAD]);
    monitor->print ("Method","sum v*v*d/t = %lg",g[INDEX_TURBULENCE_VVDoT]);
    monitor->print ("Method","sum v*v/t   = %lg",g[INDEX_TURBULENCE_VVoT]);
    monitor->print ("Method","sum v*v*d   = %lg",g[INDEX_TURBULENCE_VVD]);
    monitor->print ("Method","sum v*v     = %lg",g[INDEX_TURBULENCE_VV]);
    monitor->print ("Method","sum d*d     = %lg",g[INDEX_TURBULENCE_DD]);

    monitor->print ("Method","sum d*ax    = %lg",g[INDEX_TURBULENCE_DAx]);
    monitor->print ("Method","sum d*ay    = %lg",g[INDEX_TURBULENCE_DAy]);
    monitor->print ("Method","sum d*az    = %lg",g[INDEX_TURBULENCE_DAz]);

    monitor->print ("Method","sum d*vx    = %lg",g[INDEX_TURBULENCE_DVx]);
    monitor->print ("Method","sum d*vy    = %lg",g[INDEX_TURBULENCE_DVy]);
    monitor->print ("Method","sum d*vz    = %lg",g[INDEX_TURBULENCE_DVz]);

    monitor->print ("Method","sum d*ln(d) = %lg",g[INDEX_TURBULENCE_DlnD]);
    
    monitor->print ("Method","min d       = %lg",g[INDEX_TURBULENCE_minD]);
    monitor->print ("Method","max d       = %lg",g[INDEX_TURBULENCE_maxD]);
    monitor->print ("Method","norm        = %lg",norm);

    monitor->print ("Method","kinetic energy          %lg",
		    0.50*g[INDEX_TURBULENCE_VVD]/n);
    monitor->print ("Method","mass weighted rms Mach  %lg",
		    sqrt(g[INDEX_TURBULENCE_VVDoT]/n));
    monitor->print ("Method","volume weighed rms Mach %lg",
		    sqrt(g[INDEX_TURBULENCE_VVoT]/n));
    monitor->print ("Method","rms Velocity            %lg",
		    sqrt(g[INDEX_TURBULENCE_VV]/n));
    monitor->print ("Method","Density variance        %lg",
		    sqrt(g[INDEX_TURBULENCE_DD]/n));
    monitor->print ("Method","min/max Density         %lg",
		    g[INDEX_TURBULENCE_minD] /
		    g[INDEX_TURBULENCE_maxD]);                  
  }

  if (!comm_block->is_leaf()) return;

  const int p = field_precision (0);

  if      (p == precision_single)    
    compute_resume_<float> (comm_block,msg);
  else if (p == precision_double)    
    compute_resume_<double> (comm_block,msg);
  else if (p == precision_quadruple) 
    compute_resume_<long double> (comm_block,msg);
}

//----------------------------------------------------------------------

template <class T>
void EnzoMethodTurbulence::compute_resume_ 
(CommBlock * comm_block,
 CkReductionMsg * msg) throw()
{

  // Compute normalization

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);

  int ndx,ndy,ndz;
  domain_size(&ndx,&ndy,&ndz);
  int nd = ndx*ndy*ndz;

  double * g = enzo_block->method_turbulence_data;

  double dt = comm_block->dt();

  double norm = (edot_ != 0.0) ?
    ( sqrt(g[0]*g[0] + 2.0*nd*g[1]*dt*edot_) - g[0] ) / g[1] : 0.0;

  // ASSUMES CONSTANT TIME STEP

  double dt0 = dt;
  norm = (dt/dt0)*norm;

  const int rank = this->rank();

  Field field = enzo_block->block()->field();

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
  field.ghosts(0,&gx,&gy,&gz);
  int nbx = nx + 2*gx;
  int nby = ny + 2*gy;

  // compute bulk momentum
  const T bm[3] = 
    { T(g[INDEX_TURBULENCE_DAx]/nd),
      T(g[INDEX_TURBULENCE_DAy]/nd),
      T(g[INDEX_TURBULENCE_DAz]/nd)};

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

  enzo_block->compute_stop();
  
}
