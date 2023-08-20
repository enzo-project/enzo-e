// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialCosmology.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-07-24
/// @brief    [\ref Enzo] Declaration of the EnzoInitialCosmology class


#include "enzo.hpp"
#include "charm_enzo.hpp"

//----------------------------------------------------------------------

EnzoInitialCosmology::EnzoInitialCosmology
(int cycle, double time,
 double gamma,
 double temperature) throw()
  : Initial (cycle,time),
    gamma_(gamma),
    temperature_(temperature)
{
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();
  EnzoUnits * units = enzo::units();

  // set units to cosmology
  units->set_cosmology(cosmology);

  // set initial redshift
  double r0 = cosmology->initial_redshift();
  cosmology->set_current_redshift(r0);

  // set initial time based on redshift
  enzo::simulation()->set_time(cosmology->time_from_redshift(r0));
}

//----------------------------------------------------------------------

void EnzoInitialCosmology::enforce_block
( Block * block, const Hierarchy * hierarchy ) throw()
{

  // "If temperature is left unset, set it assuming that T=550 K at z=200"
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (temperature_ == 0.0) {
    temperature_ = 550.0 *
      pow((1.0 + cosmology->initial_redshift())/(1.0 + 200.00), 2.0);
  }

  EnzoInitialCosmology::init_cosmology(block, temperature_, gamma_);

  const GrackleChemistryData * grackle_chem = enzo::grackle_chemistry();
  
  // initialize chemistry fields if doing multispecies
  if (grackle_chem && ((grackle_chem->get<int>("metal_cooling") == 1) ||
                       (grackle_chem->get<int>("primordial_chemistry") > 0))) {

#ifdef CONFIG_USE_GRACKLE
    const EnzoMethodGrackle * grackle_method = enzo::grackle_method();
    grackle_field_data grackle_fields_; 

    //create data struct to be fed into grackle
    grackle_method->setup_grackle_fields(block, & grackle_fields_);

    //initialize density fields for various chemical species
    grackle_method->update_grackle_density_fields(block, & grackle_fields_);
#endif
  }

  block->initial_done();  
}

//----------------------------------------------------------------------

void EnzoInitialCosmology::init_cosmology
(Block * block, double temperature, double gamma) throw()
{
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  // Set initial time based on initial redshift
  block->set_time(enzo::simulation()->time());

  if (temperature < 0){ return; }

  EnzoUnits * units = enzo::units();

  const double default_mu = 0.6;

  const double internal_energy = temperature/units->kelvin_per_energy_units()
    /default_mu/(gamma-1.0);

  Field field = block->data()->field();

  // initialize energy fields
  enzo_float * ei = (enzo_float *) field.values("internal_energy");
  enzo_float * et = (enzo_float *) field.values("total_energy");
  enzo_float * vx = (enzo_float *) field.values("velocity_x");
  enzo_float * vy = (enzo_float *) field.values("velocity_y");
  enzo_float * vz = (enzo_float *) field.values("velocity_z");

  int mx,my,mz;
  field.dimensions (0,&mx,&my,&mz);

  for (int iz = 0; iz < mz; iz++) {
    for (int iy = 0; iy < my; iy++) {
      for (int ix = 0; ix < mx; ix++) {
        int i = ix + mx*(iy + my*iz);
        ei[i] = internal_energy;
        et[i] = ei[i] + 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
      }
    }
  }

}
