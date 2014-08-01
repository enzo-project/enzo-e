// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulence.cpp
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
 double mach_number)
  : Method(),
    density_initial_(density_initial),
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

  FieldBlock * field_block = comm_block->block()->field_block();

  const EnzoConfig * enzo_config = static_cast<const EnzoConfig*>
    (enzo_block->simulation()->config());

  EnzoMethodTemperature method_temperature
    (enzo_config->ppm_density_floor,
     enzo_config->ppm_temperature_floor,
     enzo_config->ppm_mol_weight);

  method_temperature.compute(enzo_block);

  enzo_float *  density = (enzo_float *) field_block->values("density");
  enzo_float *  velocity[3] = {
    (enzo_float *) field_block->values("velocity_x"),
    (enzo_float *) field_block->values("velocity_y"),
    (enzo_float *) field_block->values("velocity_z") };
  enzo_float * driving[3] = {
    (enzo_float *) field_block->values("driving_x"),
    (enzo_float *) field_block->values("driving_y"),
    (enzo_float *) field_block->values("driving_z") };
  enzo_float * temperature = (enzo_float *) field_block->values("temperature");

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);
  int gx,gy,gz;
  field_block->ghosts(0,&gx,&gy,&gz);
  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;

  const int n = sizeof(enzo_block->method_turbulence_data)/sizeof(double);
  double * g = enzo_block->method_turbulence_data;
  
  ASSERT1 ("EnzoMethodTurbulence::compute()",
	   "Size of EnzoBlock::method_turbulence_data array %d must be at least 9",
	   n, n >= 9);

  for (int i=0; i<7; i++) g[i] = 0.0;
  g[7] = std::numeric_limits<double>::max();
  g[8] = std::numeric_limits<double>::min();

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

	    g[0] += v*a*d;
	    g[1] += a*a*d;
	    g[2] += v2*d*ti;
	    g[3] += v2*ti;
	    g[4] += v2*d;
	    g[5] += v2;
	  }
	  g[6] += d*d;
	  g[7] = std::min(g[7],(double) d);
	  g[8] = std::max(g[8],(double) d);

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
  const int n = 9 * sizeof(double);
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
  for (int i=0; i<9; i++) {
    g[i] = ((double *)msg->getData())[i];
  }

  // compute normalization

  int nx,ny,nz;
  domain_size(&nx,&ny,&nz);
  int n = nx*ny*nz;

  double dt = comm_block->dt();

  // Compute edot_ if needed
  const int rank = this->rank();


  double xdm,ydm,zdm;
  lower_domain(&xdm,&ydm,&zdm);
  double xdp,ydp,zdp;
  upper_domain(&xdp,&ydp,&zdp);

  if (edot_ < 0.0) {
    double domain_x = (rank >= 1) ? (xdp - xdm) : 1.0;
    double domain_y = (rank >= 2) ? (ydp - ydm) : 1.0;
    double domain_z = (rank >= 3) ? (zdp - zdm) : 1.0;
    double box_size = domain_x;
    double box_mass = domain_x * domain_y * domain_z * density_initial_;

    float v_rms = mach_number_ / sqrt(1); //Sound speed is one.

    edot_ = 0.81/box_size*box_mass*v_rms*v_rms*v_rms;
 
    /* Approximate correction to the MacLow's factor (see eqs (7) - (8))
       for **this PPM implementation**. Seems to be OK for 64^3, 128^3 and 256^3
       Mach=3,6,10 simulations of **solenoidally** driven turbulence. */
 
    edot_  *= 0.8;
 
  }

  double norm = (edot_ != 0.0) ?
    ( sqrt(g[0]*g[0] + 2.0*n*g[1]*dt*edot_) - g[0] ) / g[1] : 0.0;

  // ASSUMES CONSTANT TIME STEP

  double dt0 = dt;
  norm = (dt/dt0)*norm;


  if (comm_block->index().is_root()) {

    Monitor * monitor = comm_block->simulation()->monitor();

    monitor->print ("Method","sum v*a*d   = %lg",g[0]);
    monitor->print ("Method","sum a*a*d   = %lg",g[1]);
    monitor->print ("Method","sum v*v*d/t = %lg",g[2]);
    monitor->print ("Method","sum v*v/t   = %lg",g[3]);
    monitor->print ("Method","sum v*v*d   = %lg",g[4]);
    monitor->print ("Method","sum v*v     = %lg",g[5]);
    monitor->print ("Method","sum d*d     = %lg",g[6]);
    monitor->print ("Method","min d       = %lg",g[7]);
    monitor->print ("Method","max d       = %lg",g[8]);
    monitor->print ("Method","norm        = %lg",norm);

    monitor->print ("Method","kinetic energy          %lg", 0.50*g[4]/n);
    monitor->print ("Method","mass weighted rms Mach  %lg",sqrt(g[2]/n));
    monitor->print ("Method","volume weighed rms Mach %lg",sqrt(g[3]/n));
    monitor->print ("Method","rms Velocity            %lg",sqrt(g[5]/n));
    monitor->print ("Method","Density variance        %lg", sqrt(g[6]/n));
    monitor->print ("Method","min/max Density         %lg",g[7]/g[8]);                  
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

  FieldBlock * field_block = enzo_block->block()->field_block();

  T * te = (T*) field_block->values ("total_energy");
  T * v3[3] = { (T*) field_block->values ("velocity_x"),
		(T*) field_block->values ("velocity_y"),
		(T*) field_block->values ("velocity_z") };
  T * a3[3] = { (T*) field_block->values ("driving_x"),
		(T*) field_block->values ("driving_y"),
		(T*) field_block->values ("driving_z") };

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);
  int gx,gy,gz;
  field_block->ghosts(0,&gx,&gy,&gz);
  int nbx = nx + 2*gx;
  int nby = ny + 2*gy;

  for (int iz=gz; iz<gz+nz; iz++) {
    for (int iy=gy; iy<gy+ny; iy++) {
      for (int ix=gx; ix<gx+nx; ix++) {
	int i = ix + nbx*(iy + nby*iz);
	for (int id=0; id<rank; id++) {
	  te[i] += v3[id][i]*a3[id][i]*norm + 0.5*a3[id][i]*norm*a3[id][i]*norm;
	  v3[id][i] += a3[id][i]*norm;
	}
      }
    }
  }
}
