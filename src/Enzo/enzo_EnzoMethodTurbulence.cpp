// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTurbulence.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Wed Jul 23 00:31:04 UTC 2014
/// @brief    Implements the EnzoMethodTurbulence class

#include "cello.hpp"

#include "enzo.hpp"

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
  int ndz = nz + 2*gz;

  const int n = sizeof(enzo_block->method_turbulence_data)/sizeof(long double);
  long double * g = enzo_block->method_turbulence_data;
  
  ASSERT1 ("EnzoMethodTurbulence::compute()",
	   "Size of EnzoBlock::method_turbulence_data array %d must be at least 9",
	   n, n >= 9);

  g[7] = std::numeric_limits<double>::min();
  g[8] = std::numeric_limits<double>::max();

  const int rank = comm_block->rank();

  for (int id=0; id<rank; id++) {

    for (int iz=gz; iz<gz+nz; iz++) {
      for (int iy=gy; iy<gy+ny; iy++) {
	for (int ix=gx; ix<gx+nx; ix++) {

	  int i = ix + ndx*(iy + ndy*iz);

	  enzo_float d  = density[i];
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
	  g[6] += d*d;
	  g[7] = std::min(g[7],(long double) d);
	  g[8] = std::max(g[8],(long double) d);

	}
      }
    }
  }

  enzo_block->method_turbulence_begin();
}

//----------------------------------------------------------------------

void EnzoBlock::method_turbulence_begin()
{
  const int n = 7 * sizeof(long double);
  long double * g = method_turbulence_data;
  //  CkCallback cb(CkIndex_EnzoBlock::p_method_turbulence_sum(NULL),  thisProxy);
  //  contribute(n,&g[0],CkReduction::sum_long_double,cb);
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_turbulence_sum(CkReductionMsg *msg)
{
  
  const int n = 1 * sizeof(double);
  long double * g = method_turbulence_data;
  //  CkCallback cb(CkIndex_EnzoBlock::p_method_turbulence_min(NULL),  thisProxy);
  //  contribute(n,&g[7],CkReduction::sum_long_double,cb);
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_turbulence_min(CkReductionMsg *msg)
{
  const int n = 1 * sizeof(double);
  long double * g = method_turbulence_data;
  //  CkCallback cb(CkIndex_EnzoBlock::p_method_turbulence_max(NULL),  thisProxy);
  //  contribute(n,&g[8],CkReduction::min_long_double,cb);
}

//----------------------------------------------------------------------

void EnzoBlock::p_method_turbulence_max(CkReductionMsg *msg)
{
  //
}

//----------------------------------------------------------------------

void EnzoBlock::method_turbulence_end()
{
}
