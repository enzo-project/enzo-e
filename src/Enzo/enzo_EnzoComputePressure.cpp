// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputePressure.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoComputePressure class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputePressure::EnzoComputePressure (double gamma,
					  int comoving_coordinates)
  : Compute(),
    gamma_(gamma),
    comoving_coordinates_(comoving_coordinates)
{
}

//----------------------------------------------------------------------

void EnzoComputePressure::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);

  p | gamma_;
  p | comoving_coordinates_;

}

//----------------------------------------------------------------------

void EnzoComputePressure::compute ( Block * block) throw()
{
  Field field = block->data()->field();
  if (field.precision(0) == precision_single) {
    compute_<float>(block);
  } else if (field.precision(0) == precision_double) {
    compute_<double>(block);
  }
}

//----------------------------------------------------------------------

template <typename T>
void EnzoComputePressure::compute_(Block * block)
{

  if (!block->is_leaf()) return;

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = enzo_block->data()->field();

  T * p = (T*) field.values("pressure");
  T * d = (T*) field.values("density");
  const int rank = enzo_block->rank();
  T * v3[3] = 
    { (T*) ((rank >= 1) ? field.values("velocity_x") : NULL),
      (T*) ((rank >= 2) ? field.values("velocity_y") : NULL),
      (T*) ((rank >= 3) ? field.values("velocity_z") : NULL) };

  T * te = (T*) field.values("total_energy");

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  int gx,gy,gz;
  field.ghosts (0,&gx,&gy,&gz);
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

