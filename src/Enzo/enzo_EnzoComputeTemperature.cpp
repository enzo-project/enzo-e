// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeTemperature.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2014-10-27 22:37:41
/// @brief    Implements the EnzoComputeTemperature class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoComputeTemperature::EnzoComputeTemperature 
(double density_floor,
 double temperature_floor,
 double mol_weight,
 int comoving_coordinates ) 
  : Compute(),
    density_floor_(density_floor),
    temperature_floor_(temperature_floor),
    mol_weight_(mol_weight),
    comoving_coordinates_(comoving_coordinates)
{
}

//----------------------------------------------------------------------

void EnzoComputeTemperature::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Compute::pup(p);

  p | density_floor_;
  p | temperature_floor_;
  p | mol_weight_;
  p | comoving_coordinates_;

}

//----------------------------------------------------------------------

void EnzoComputeTemperature::compute ( Block * block) throw()
{

  if (!block->is_leaf()) return;

  Field field = block->data()->field();

  if (field.precision(0) == precision_single) {
    compute_<float>(block);
  } else if (field.precision(0) == precision_double) {
    compute_<double>(block);
  }
}

//----------------------------------------------------------------------

template <typename T>
void EnzoComputeTemperature::compute_(Block * block)
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = enzo_block->data()->field();

  EnzoComputePressure compute_pressure(EnzoBlock::Gamma,
				       comoving_coordinates_);

  compute_pressure.compute(block);

  T * t = (T*) field.values("temperature");
  T * d = (T*) field.values("density");
  T * p = (T*) field.values("pressure");

  const int rank = block->rank();

  int nx,ny,nz;
  field.size(&nx,&ny,&nz);

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);
  if (rank < 1) gx = 0;
  if (rank < 2) gy = 0;
  if (rank < 3) gz = 0;

  int m = (nx+2*gx) * (ny+2*gy) * (nz+2*gz);

  for (int i=0; i<m; i++) {
    T density     = std::max(d[i], (T) density_floor_);
    T temperature = p[i] * mol_weight_ / density;
    t[i] = std::max(temperature, (T)temperature_floor_);
  }
}

