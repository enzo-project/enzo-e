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

#ifndef CONFIG_USE_GRACKLE

  for (int i=0; i<m; i++) {
    enzo_float density     = std::max(d[i], (enzo_float) density_floor_);
    enzo_float temperature = p[i] * mol_weight_ / density;
    t[i] = std::max(temperature, (enzo_float)temperature_floor_);
  }

#else

  if (grackle_data->primordial_chemistry == 0){

    for (int i=0; i<m; i++) {
      enzo_float density     = std::max(d[i], (enzo_float) density_floor_);
      enzo_float temperature = p[i] * mol_weight_ / density;
      t[i] = std::max(temperature, (enzo_float)temperature_floor_);
    }

  } else{
    enzo_float inv_metal_mol = 1.0 / cello::mu_metal;

    enzo_float * HI_density    = (enzo_float*) field.values("HI_density");
    enzo_float * HII_density   = (enzo_float*) field.values("HII_density");
    enzo_float * HeI_density   = (enzo_float*) field.values("HeI_density");
    enzo_float * HeII_density  = (enzo_float*) field.values("HeII_density");
    enzo_float * HeIII_density = (enzo_float*) field.values("HeIII_density");
    enzo_float * H2I_density   = field.is_field("H2I_density") ?
                                   (enzo_float*) field.values("H2I_density") : NULL;
    enzo_float * H2II_density  = field.is_field("H2II_density") ?
                                   (enzo_float*) field.values("H2II_density") : NULL;
    enzo_float * HM_density    = field.is_field("HM_density") ?
                                   (enzo_float*) field.values("HM_density") : NULL;
    enzo_float * e_density     = (enzo_float*) field.values("e_density");
    enzo_float * metal_density = field.is_field("metal_density") ?
                             (enzo_float*) field.values("metal_density") : NULL;

    for (int i=0; i<m; i++){
      enzo_float number_density
        = 0.25*(HeI_density[i] + HeII_density[i] + HeIII_density[i]) +
           HI_density[i] + HII_density[i] + e_density[i];

      if(grackle_data->primordial_chemistry > 1){
        number_density += HM_density[i] +
          0.5 * (H2I_density[i] + H2II_density[i]);
      }

      // if metals exist
       if( metal_density != NULL)
         number_density += metal_density[i] * inv_metal_mol;

       // use an approximate density floor for number density
       number_density = std::max(number_density,
                        enzo_float(density_floor_ / mol_weight_));
       enzo_float temperature = p[i] / number_density;
       t[i] = std::max(temperature, (enzo_float) temperature_floor_);
    }


  } // end if primordial_chemistry

#endif // CONFIG_USE_GRACKLE
}
