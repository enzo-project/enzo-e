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
  EnzoUnits * units = enzo::units();

  // "If temperature is left unset, set it assuming that T=550 K at z=200"
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  if (temperature_ == 0.0) {
    temperature_ = 550.0 *
      pow((1.0 + cosmology->initial_redshift())/(1.0 + 200.00), 2.0);
  }
  
  block->set_time(enzo::simulation()->time());
  
  const double default_mu = 0.6;

  const double internal_energy = temperature_/units->kelvin_per_energy_units()
    /default_mu/(gamma_-1.0);

  Field field = block->data()->field();
  
  enzo_float * ei = (enzo_float *) field.values("internal_energy");
  enzo_float * et = (enzo_float *) field.values("total_energy");
  enzo_float * vx = (enzo_float *) field.values("velocity_x");
  enzo_float * vy = (enzo_float *) field.values("velocity_y");
  enzo_float * vz = (enzo_float *) field.values("velocity_z");

  int mx,my,mz;
  int gx,gy,gz;
  field.dimensions (0,&mx,&my,&mz);
  field.ghost_depth(0,&gx,&gy,&gz);

  gx=gy=gz=0;
  for (int iz=gz; iz<mz-gz; iz++) {
    for (int iy=gy; iy<my-gy; iy++) {
      for (int ix=gx; ix<mx-gx; ix++) {
	int i = ix + mx*(iy + my*iz);
	ei[i] = internal_energy;
	et[i] = ei[i] + 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
      }
    }
  }

  block->initial_done();

#ifdef CONFIG_USE_GRACKLE
  // initialize chemistry fields if doing multispecies 
  if (enzo::config()->method_grackle_chemistry)
  {
    chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

    const EnzoMethodGrackle * grackle_method = enzo::grackle_method();
 
    if ( (grackle_chemistry->primordial_chemistry > 0) || (grackle_chemistry->metal_cooling == 1))
    {
        EnzoBlock * enzo_block = enzo::block(block);
        grackle_field_data grackle_fields_; 
 
        //create data struct to be fed into grackle
        grackle_method->setup_grackle_fields(enzo_block, & grackle_fields_);

        //initialize density fields for various chemical species
        grackle_method->update_grackle_density_fields(enzo_block, & grackle_fields_);
    }
  }
#endif 

}
