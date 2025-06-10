// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodComovingExpansion.cpp
/// @author   Britton Smith (bds006@ucsd.edu)
/// @date     Wed May 24 12:25:56 PDT 2017
/// @brief    Implements comoving expansion class

#include "Enzo/cosmology/cosmology.hpp"
#include "Enzo/enzo.hpp"
//----------------------------------------------------------------------

EnzoMethodComovingExpansion::EnzoMethodComovingExpansion
( bool comoving_coordinates )
  : Method(),
    comoving_coordinates_(comoving_coordinates)
{
  cello::simulation()->refresh_set_name(ir_post_,name());

  const int rank = cello::rank();

  cello::define_field ("density");
  cello::define_field ("total_energy");
  cello::define_field ("internal_energy");
  cello::define_field ("pressure");
  if (rank >= 1) cello::define_field ("velocity_x");
  if (rank >= 2) cello::define_field ("velocity_y");
  if (rank >= 3) cello::define_field ("velocity_z");

  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_field("density");
  refresh->add_field("total_energy");
  refresh->add_field("internal_energy");
  refresh->add_field("velocity_x");
  refresh->add_field("velocity_y");
  refresh->add_field("velocity_z");

  if ( ! comoving_coordinates_ ) {
    WARNING
      ("EnzoMethodComovingExpansion::EnzoMethodComovingExpansion()",
       "Including \"comoving_expansion\" method but cosmology is disabled");
  }
}

//----------------------------------------------------------------------

void EnzoMethodComovingExpansion::compute ( Block * block) throw()
{

  if (block->state()->cycle() == 0) {
    // skip first cycle
    block->compute_done();
    return;
  }
  EnzoBlock * enzo_block = enzo::block(block);
  Field field = enzo_block->data()->field();

  Monitor * monitor = cello::monitor();
  if (block->index().is_root()) {
    monitor->print("Method", "%s redshift %.8f",
		   this->name().c_str(),
		   enzo::cosmology()->current_redshift());
  }


  /* Only do this if
     1. this is a leaf block
     2. we are using comoving coordinates
     3. baryon fields are present.
  */

  if ((block->is_leaf() &&
       comoving_coordinates_ &&
       field.field_count() > 0))
    {

      EnzoPhysicsCosmology * cosmology = enzo::cosmology();

      ASSERT ("EnzoMethodComovingExpansion::compute()",
  	      "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
  	      ! (comoving_coordinates_ && (cosmology == NULL)) );

      /* Compute adot/a at time = t-1/2dt (time-centered). */

      int has_history = ((field.num_history() > 0) &&
  			 (field.history_time(1) > 0.));
      enzo_float compute_time;
      if (has_history) {
  	compute_time = 0.5 * (enzo_block->state()->time() +
  			      field.history_time(1));
      }
      else {
  	compute_time = enzo_block->state()->time();
      }

      //      printf ("DEBUG_VELOCITY time old new = %g %g\n",field.history_time(1),enzo_block->state()->time());
      enzo_float cosmo_a=1.0;
      enzo_float cosmo_dadt=0.0;
      cosmology->compute_expansion_factor (&cosmo_a, &cosmo_dadt, compute_time);
      double dt = block->state()->dt();
      enzo_float Coefficient = dt*cosmo_dadt/cosmo_a;


      /* Determine the size of the block. */

      int mx, my, mz, m, rank;
      field.dimensions(0,&mx,&my,&mz);
      m = mx*my*mz;
      rank = cello::rank();

      /* If we can, compute the pressure at the mid-point.
      	 We can, because we will always have an old baryon field now. */
      const int in = cello::index_static();
      enzo_float gamma = enzo::fluid_props()->gamma();

      // hard-code hydromethod for PPM for now
      int HydroMethod = 0;

      // hard-code CR method to off
      int CRModel = 0;
      enzo_float * cr_field_new = NULL;
      enzo_float * cr_field_old = NULL;

      /* Get the necessary fields.
  	 field.values(<field_name>, 0) is the field at the current time.
  	 field.values(<field_name>, 1) is the field at the previous time.
      */

      int i_new = 0;
      int i_old = has_history ? 1 : 0;

      enzo_float * density_new         =
  	(enzo_float *) field.values("density", i_new);
      enzo_float * density_old         =
  	(enzo_float *) field.values("density", i_old);
      enzo_float * total_energy_new    =
  	(enzo_float *) field.values("total_energy", i_new);
      enzo_float * total_energy_old    =
  	(enzo_float *) field.values("total_energy", i_old);
      enzo_float * internal_energy_new =
  	(enzo_float *) field.values("internal_energy", i_new);
      enzo_float * internal_energy_old =
  	(enzo_float *) field.values("internal_energy", i_old);

      enzo_float * velocity_x_new = NULL;
      enzo_float * velocity_y_new = NULL;
      enzo_float * velocity_z_new = NULL;
      enzo_float * velocity_x_old = NULL;
      enzo_float * velocity_y_old = NULL;
      enzo_float * velocity_z_old = NULL;

      velocity_x_new   = (enzo_float *) field.values("velocity_x", i_new);
      velocity_x_old   = (enzo_float *) field.values("velocity_x", i_old);
      if (rank >= 2) {
  	velocity_y_new = (enzo_float *) field.values("velocity_y", i_new);
  	velocity_y_old = (enzo_float *) field.values("velocity_y", i_old);
      }
      if (rank >= 3) {
  	velocity_z_new = (enzo_float *) field.values("velocity_z", i_new);
  	velocity_z_old = (enzo_float *) field.values("velocity_z", i_old);
      }

      // Compute the pressure *now*
      enzo_float * pressure_now     = (enzo_float *) field.values("pressure");
      EnzoComputePressure compute_pressure (gamma, comoving_coordinates_);
      compute_pressure.compute(block, pressure_now);

      // If history is present, compute time-centered pressure
      enzo_float * pressure = NULL;

      if (has_history) {
        EnzoComputePressure compute_pressure_old (gamma, comoving_coordinates_);
        compute_pressure_old.set_history(i_old);

        pressure = new enzo_float[m];

        compute_pressure_old.compute(block, pressure);

        // now compute the time-centered average of the two
        for (int i = 0; i < m; i ++){
          pressure[i] = 0.5*(pressure_now[i] + pressure[i]);
        }

      } else {

        // if not, just use current pressure
        pressure = pressure_now;
      }

      int idual = static_cast<int>
        (! enzo::fluid_props()->dual_energy_config().is_disabled());
      /* Call fortran routine to do the real work. */

      FORTRAN_NAME(expand_terms)
	(
	 &rank, &m, &idual, &Coefficient,
	 (int*) &HydroMethod, &gamma,
	 pressure,
	 density_new, total_energy_new, internal_energy_new,
	 velocity_x_new, velocity_y_new, velocity_z_new,
	 density_old, total_energy_old, internal_energy_old,
	 velocity_x_old, velocity_y_old, velocity_z_old,
	 &CRModel, cr_field_new, cr_field_old);


      if (has_history){
	delete [] pressure;
	pressure = NULL;
      }
    }

  block->compute_done();
}

//----------------------------------------------------------------------

double EnzoMethodComovingExpansion::timestep( Block * block ) throw()
{

  enzo_float dtExpansion = ENZO_HUGE_VAL;

  if (!comoving_coordinates_)
    return (double) dtExpansion;

  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  ASSERT ("EnzoMethodComovingExpansion::timestep()",
	  "comoving_coordinates enabled but missing EnzoPhysicsCosmology",
	  ! (comoving_coordinates_ && (cosmology == NULL)) );

  EnzoBlock * enzo_block = enzo::block(block);

  cosmology->compute_expansion_timestep(&dtExpansion,
                                        (enzo_float) enzo_block->state()->time());

  return (double) dtExpansion;

}

//----------------------------------------------------------------------

void EnzoMethodComovingExpansion::pup (PUP::er &p)
{
  TRACEPUP;
  // NOTE: change this function whenever attributes change

  Method::pup(p);

  p | comoving_coordinates_;
}

//======================================================================
