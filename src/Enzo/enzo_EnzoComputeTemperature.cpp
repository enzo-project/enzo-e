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
 bool comoving_coordinates ) 
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

  compute_(block);
}

//----------------------------------------------------------------------

void EnzoComputeTemperature::compute_(Block * block)
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  Field field = enzo_block->data()->field();

  const int in = cello::index_static();

  EnzoComputePressure compute_pressure(EnzoBlock::Gamma[in],
				       comoving_coordinates_);

  compute_pressure.compute(block);

  enzo_float * t = (enzo_float*) field.values("temperature");
  enzo_float * d = (enzo_float*) field.values("density");
  enzo_float * p = (enzo_float*) field.values("pressure");

  int mx,my,mz;
  field.dimensions(0,&mx,&my,&mz);

  const int m = mx*my*mz;
  
  for (int i=0; i<m; i++) {
    enzo_float density     = std::max(d[i], (enzo_float) density_floor_);
    enzo_float temperature = p[i] * mol_weight_ / density;
    t[i] = std::max(temperature, (enzo_float)temperature_floor_);
  }
}

