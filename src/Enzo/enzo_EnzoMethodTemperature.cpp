// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodTemperature.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodTemperature class

#include "cello.hpp"

#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodTemperature::EnzoMethodTemperature 
(double density_floor,
 double temperature_floor,
 double mol_weight) :
  Method(),
  density_floor_(density_floor),
  temperature_floor_(temperature_floor),
  mol_weight_(mol_weight)
{
}

//----------------------------------------------------------------------

void EnzoMethodTemperature::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | density_floor_;
  p | temperature_floor_;
  p | mol_weight_;

}

//----------------------------------------------------------------------

void EnzoMethodTemperature::compute ( CommBlock * comm_block) throw()
{

  if (!comm_block->is_leaf()) return;

  initialize_(comm_block);

  if (field_precision(0) == precision_single) {
    compute_<float>(comm_block);
  } else if (field_precision(0) == precision_double) {
    compute_<double>(comm_block);
  }
}

//----------------------------------------------------------------------

template <typename T>
void EnzoMethodTemperature::compute_(CommBlock * comm_block)
{
  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (comm_block);

  FieldBlock * field_block = enzo_block->block()->field_block();

  EnzoMethodPressure method_pressure(EnzoBlock::Gamma);

  method_pressure.compute(comm_block);

  T * t = (T*) field_block->values("temperature");
  T * d = (T*) field_block->values("density");
  T * p = (T*) field_block->values("pressure");

  const int rank = this->rank();

  int nx,ny,nz;
  field_block->size(&nx,&ny,&nz);

  int gx,gy,gz;
  field_block->ghosts (0,&gx,&gy,&gz);
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

