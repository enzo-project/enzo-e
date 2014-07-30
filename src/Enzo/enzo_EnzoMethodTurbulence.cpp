// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulence.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:31:04 UTC 2014
/// @brief    Implements the EnzoMethodTurbulence class

#include "cello.hpp"

#include "enzo.hpp"

#include "enzo.decl.h"

//----------------------------------------------------------------------

EnzoMethodTurbulence::EnzoMethodTurbulence () : Method()
{
  // TURBULENCE parameters initialized in EnzoBlock::initialize()
}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodTurbulence::compute ( CommBlock * comm_block) throw()
{
  TRACE("EnzoMethodTurbulence::compute()");

  //  INCOMPLETE("EnzoMethodTurbulence::compute()");

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);

  FieldBlock * field_block = comm_block->block()->field_block();

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

  const int rank = comm_block->rank();

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

  enzo_block->method_turbulence_begin();
}

//----------------------------------------------------------------------

extern CkReduction::reducerType r_method_turbulence_type;

void EnzoBlock::method_turbulence_begin()
{
  TRACE("EnzoMethodTurbulence::r_method_turbulence_type()");
  const int n = 9 * sizeof(double);
  double * g = method_turbulence_data;
  CkCallback cb (CkIndex_EnzoBlock::p_method_turbulence_end(NULL),thisProxy);
  contribute(n,&g,r_method_turbulence_type,cb);
}

//----------------------------------------------------------------------

// SEE main.cpp for implementation

// CkReductionMsg * r_method_turbulence(int n, CkReductionMsg ** msgs)

//----------------------------------------------------------------------

void EnzoBlock::p_method_turbulence_end(CkReductionMsg * msg)
{
  TRACE("EnzoMethodTurbulence::p_method_turbulence_end()");
  double * g = method_turbulence_data;
  if (((CommBlock *)this)->index().is_root()) {
    simulation()->monitor()->print ("Method","sum v*a*d   = %lg",g[0]);
    simulation()->monitor()->print ("Method","sum a*a*d   = %lg",g[1]);
    simulation()->monitor()->print ("Method","sum v*v*d/t = %lg",g[2]);
    simulation()->monitor()->print ("Method","sum v*v/t   = %lg",g[3]);
    simulation()->monitor()->print ("Method","sum v*v*d   = %lg",g[4]);
    simulation()->monitor()->print ("Method","sum v*v     = %lg",g[5]);
    simulation()->monitor()->print ("Method","sum d*d     = %lg",g[6]);
    simulation()->monitor()->print ("Method","min d       = %lg",g[7]);
    simulation()->monitor()->print ("Method","max d       = %lg",g[8]);

    int nx,ny,nz;
    simulation()->hierarchy()->root_size(&nx,&ny,&nz);
    int n = nx*ny*nz;

    simulation()->monitor()->print 
      ("Method","kinetic energy          %lg", 0.50*g[4]/n);

    simulation()->monitor()->print 
      ("Method","mass weighted rms Mach  %lg",sqrt(g[2]/n));
    simulation()->monitor()->print
      ("Method","volume weighed rms Mach %lg",sqrt(g[3]/n));
    simulation()->monitor()->print
      ("Method","rms Velocity            %lg",sqrt(g[5]/n));
    simulation()->monitor()->print
      ("Method","Density variance        %lg", sqrt(g[6]/n));
    simulation()->monitor()->print
      ("Method","min/max Density         %lg",g[7]/g[8]);                  
  }
}
