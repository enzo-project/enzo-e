// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeTemperature.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
///           Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-05-07
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

  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  enzo_float * t = field.is_field("temperature") ?
                    (enzo_float*) field.values("temperature", i_hist_) : NULL;

  if (!t) {
    ERROR("EnzoComputeTemperature::compute()",
          " 'temperature' field is not defined as a permanent field");
  }

  compute(block, t);
}

//---------------------------------------------------------------------

void EnzoComputeTemperature::compute ( Block * block, enzo_float * t) throw()
{

  if (!block->is_leaf()) return;

  compute_(block, t);
}

//----------------------------------------------------------------------

void EnzoComputeTemperature::compute_(Block * block,
                                      enzo_float * t,
                                      bool recompute_pressure /* true */
#ifdef CONFIG_USE_GRACKLE
                                    , code_units * grackle_units /* NULL */ ,
                                      grackle_field_data * grackle_fields /* NULL */
#endif
                                    )
{
  EnzoBlock * enzo_block = enzo::block(block);

  Field field = enzo_block->data()->field();

  const int in = cello::index_static();

  if (enzo::config()->method_grackle_use_grackle){

#ifdef CONFIG_USE_GRACKLE
    code_units grackle_units_;
    grackle_field_data grackle_fields_;

    // setup grackle units if they are not already provided
    if (!grackle_units){
      grackle_units = &grackle_units_;
      EnzoMethodGrackle::setup_grackle_units(enzo_block, grackle_units, i_hist_);
    }

    // if grackle fields are not provided, define them
    bool delete_grackle_fields = false;
    if (!grackle_fields){
      grackle_fields  = &grackle_fields_;
      EnzoMethodGrackle::setup_grackle_fields(enzo_block, grackle_fields, i_hist_);
      delete_grackle_fields = true;
    }

    // Note: compute pressure is handled internally in Grackle.
    //       Do not need to recompute this here

    // temperature is returned in units of K
    if (calculate_temperature(grackle_units, grackle_fields, t) == ENZO_FAIL){
      ERROR("EnzoComputeTemperature::compute_()",
            "Error in call to Grackle's compute_temperature routine.\n");
    }

    if (delete_grackle_fields){
      EnzoMethodGrackle::delete_grackle_fields(grackle_fields);
    }
#else
    ERROR("EnzoComputeTemperature::compute_()",
          "Attempting to compute temperature with method Grackle "
          "but Enzo-E has not been compiled with Grackle (set use_grackle = 1) \n");
#endif

  } else {

    EnzoUnits * enzo_units = enzo::units();

    int mx,my,mz;
    field.dimensions(0,&mx,&my,&mz);

    const int m = mx*my*mz;

    enzo_float * d = (enzo_float*) field.values("density", i_hist_);
    enzo_float * p = (enzo_float*) field.values("pressure", i_hist_);

    EnzoComputePressure compute_pressure(EnzoBlock::Gamma[in],
                                       comoving_coordinates_);
    compute_pressure.set_history(i_hist_);

    if (recompute_pressure) compute_pressure.compute(block, p);

    for (int i=0; i<m; i++) {
      enzo_float density     = std::max(d[i], (enzo_float) density_floor_);
      enzo_float temperature = p[i] * mol_weight_ / density;
      t[i] = std::max(temperature, (enzo_float)temperature_floor_) * enzo_units->temperature();
    }
  }

  return;
}
