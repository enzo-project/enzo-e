// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodPressure.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodPressure class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodPressure::EnzoMethodPressure (double gamma) :
  Method(),gamma_(gamma)
{
}

//----------------------------------------------------------------------

void EnzoMethodPressure::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

}

//----------------------------------------------------------------------

void EnzoMethodPressure::compute ( CommBlock * comm_block) throw()
{
  initialize_(comm_block);

  if (field_precision(0) == precision_single) {
    compute_<float>(comm_block);
  } else if (field_precision(0) == precision_double) {
    compute_<double>(comm_block);
  }
}

//----------------------------------------------------------------------

template <typename T>
void EnzoMethodPressure::compute_(CommBlock * comm_block)
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);

  FieldBlock * field_block = enzo_block->block()->field_block();

  T * p = (T*) field_block->values("pressure");
  T * d = (T*) field_block->values("density");
  T * v3[3] = 
    { (T*)field_block->values("velocity_x"),
      (T*)field_block->values("velocity_y"),
      (T*)field_block->values("velocity_z")  };

  T * te = (T*)field_block->values("total_energy");

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  int gx,gy,gz;
  field_block->ghosts (0,&gx,&gy,&gz);
  const int rank = this->rank();
  if (rank < 1) gx = 0;
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  int m = (nx+2*gx) * (ny+2*gy) * (nz+2*gz);
  T gm1 = gamma_ - 1.0;
  for (int i=0; i<m; i++) {
    T e= te[i];
    if (rank >= 1) e -= 0.5*v3[0][i]*v3[0][i];
    if (rank >= 2) e -= 0.5*v3[1][i]*v3[1][i];
    if (rank >= 3) e -= 0.5*v3[2][i]*v3[2][i];
    p[i] = gm1 * d[i] * e;
  }
}

