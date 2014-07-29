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

  enzo_float *  d = (enzo_float *) field_block->values("density");
  enzo_float *  v3[3] = {
    (enzo_float *) field_block->values("velocity_x"),
    (enzo_float *) field_block->values("velocity_y"),
    (enzo_float *) field_block->values("velocity_z") };
  enzo_float * a3[3] = {
    (enzo_float *) field_block->values("driving_x"),
    (enzo_float *) field_block->values("driving_y"),
    (enzo_float *) field_block->values("driving_z") };
  enzo_float * t = (enzo_float *) field_block->values("temperature");

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);
  int gx,gy,gz;
  field_descr->ghosts(0,&gx,&gy,&gz);

  int ndx = nx + 2*gx;
  int ndy = ny + 2*gy;
  int ndz = nz + 2*gz;

  long double g[9] = {0.0};
  g[7] = std::numeric_limits<double>::min();
  g[8] = std::numeric_limits<double>::max();

  for (int iz=gz; iz<gz+nz; iz++) {
    for (int iy=gy; iy<gy+ny; iy++) {
      for (int ix=gx; ix<gx+nx; ix++) {
	int i = ix + ndx*(iy + ndy*iz);
	g[0] += v3[id][i] * a3[id][i] * d[i];
      }
    }
  }

}

