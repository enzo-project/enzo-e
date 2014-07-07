// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodGrackle class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------

EnzoMethodGrackle::EnzoMethodGrackle 
(EnzoConfig * config)
  : Method()
#ifdef CONFIG_USE_GRACKLE
  , chemistry_(),
    units_()
#endif /* CONFIG_USE_GRACKLE */

{
#ifdef CONFIG_USE_GRACKLE

  /// Initialize parameters

  units_     = config->method_grackle_units;
  chemistry_ = config->method_grackle_chemistry;

  const gr_float a_value = 
    1. / (1. + config->physics_cosmology_initial_redshift);

  if (initialize_chemistry_data
      (chemistry_, units_, a_value) == 0) {
    ERROR("EnzoMethodGrackle::EnzoMethodGrackle()",
	  "Error in initialize_chemistry_data");
  }
  
#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

void EnzoMethodGrackle::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

  TRACEPUP;

#ifdef CONFIG_USE_GRACKLE

  Method::pup(p);

  p | chemistry_;
  p | units_;

#endif /* CONFIG_USE_GRACKLE */

}

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute ( CommBlock * comm_block) throw()
{
#ifdef CONFIG_USE_GRACKLE
  printf ("EnzoMethodGrackle::compute()\n");
  initialize_(comm_block);

  // ASSUMES ALL ARRAYS ARE THE SAME SIZE
  gr_int m[3];
  array_dimension (0,m,m+1,m+2);

  int gx,gy,gz;
  ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  array_size (&nx,&ny,&nz);

  gr_int grid_start[3] = {gx,      gy,      gz};
  gr_int grid_end[3]   = {gx+nx-1, gy+ny-1, gz+nz-1} ;

  // ASSUMES COSMOLOGY = false
  double a_value = 1.0;

  gr_int rank = this->rank();

  gr_float * density = (gr_float *) field_array(field_id("density"));
  gr_float * energy  = (gr_float *) field_array(field_id("energy"));
  gr_float * x_velocity = (gr_float *) field_array(field_id("x_velocity"));
  gr_float * y_velocity = (gr_float *) field_array(field_id("y_velocity"));
  gr_float * z_velocity = (gr_float *) field_array(field_id("z_velocity"));
  gr_float * HI_density = (gr_float *) field_array(field_id("HI_density"));
  gr_float * HII_density = (gr_float *) field_array(field_id("HII_density"));
  gr_float * HM_density = (gr_float *) field_array(field_id("HM_density"));
  gr_float * HeI_density = (gr_float *) field_array(field_id("HeI_density"));
  gr_float * HeII_density = (gr_float *) field_array(field_id("HeII_density"));
  gr_float * HeIII_density = (gr_float *) field_array(field_id("HeIII_density"));
  gr_float * H2I_density = (gr_float *) field_array(field_id("H2I_density"));
  gr_float * H2II_density = (gr_float *) field_array(field_id("H2II_density"));
  gr_float * DI_density = (gr_float *) field_array(field_id("DI_density"));
  gr_float * DII_density = (gr_float *) field_array(field_id("DII_density"));
  gr_float * HDI_density = (gr_float *) field_array(field_id("HDI_density"));
  gr_float * e_density = (gr_float *) field_array(field_id("e_density"));
  gr_float * metal_density = (gr_float *) field_array(field_id("metal_density"));

  double dt = time_step();

  printf ("density %p e_density %p\n",density,e_density);
  if (solve_chemistry
      (chemistry_, units_,
       a_value, dt,
       rank, m,
       grid_start, grid_end,
       density, 
       energy,
       x_velocity,
       y_velocity,
       z_velocity,
       HI_density,
       HII_density,
       HM_density,
       HeI_density,
       HeII_density,
       HeIII_density,
       H2I_density,
       H2II_density,
       DI_density,
       DII_density,
       HDI_density,
       e_density,
       metal_density) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in solve_chemistry");
  }

  printf ("%s:%d TRACE exited solve_chemistry()\n",
	  __FILE__,__LINE__); fflush(stdout);

  if (calculate_cooling_time
      (chemistry_, units_,
       a_value,
       rank, m,
       grid_start, grid_end,
       (gr_float *) field_array(field_id("density")), 
       (gr_float *) field_array(field_id("energy")),
       (gr_float *) field_array(field_id("x_velocity")), 
       (gr_float *) field_array(field_id("y_velocity")), 
       (gr_float *) field_array(field_id("z_velocity")),
       (gr_float *) field_array(field_id("HI_density")), 
       (gr_float *) field_array(field_id("HII_density")), 
       (gr_float *) field_array(field_id("HM_density")),
       (gr_float *) field_array(field_id("HeI_density")), 
       (gr_float *) field_array(field_id("HeII_density")), 
       (gr_float *) field_array(field_id("HeIII_density")),
       (gr_float *) field_array(field_id("H2I_density")), 
       (gr_float *) field_array(field_id("H2II_density")),
       (gr_float *) field_array(field_id("DI_density")), 
       (gr_float *) field_array(field_id("DII_density")), 
       (gr_float *) field_array(field_id("HDI_density")),
       (gr_float *) field_array(field_id("e_density")), 
       (gr_float *) field_array(field_id("metal_density")), 
       (gr_float *) field_array(field_id("cooling_time"))) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in calculate_cooling_time.\n");
  }

  if (calculate_temperature
      (chemistry_, units_,
       rank, m,
       (gr_float *) field_array(field_id("density")), 
       (gr_float *) field_array(field_id("energy")),
       (gr_float *) field_array(field_id("HI_density")), 
       (gr_float *) field_array(field_id("HII_density")), 
       (gr_float *) field_array(field_id("HM_density")),
       (gr_float *) field_array(field_id("HeI_density")), 
       (gr_float *) field_array(field_id("HeII_density")), 
       (gr_float *) field_array(field_id("HeIII_density")),
       (gr_float *) field_array(field_id("H2I_density")), 
       (gr_float *) field_array(field_id("H2II_density")),
       (gr_float *) field_array(field_id("DI_density")), 
       (gr_float *) field_array(field_id("DII_density")), 
       (gr_float *) field_array(field_id("HDI_density")),
       (gr_float *) field_array(field_id("e_density")), 
       (gr_float *) field_array(field_id("metal_density")), 
       (gr_float *) field_array(field_id("temperature"))) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in calculate_temperature.\n");
  }

  if (calculate_pressure
      (chemistry_, units_,
       rank, m,
       (gr_float *) field_array(field_id("density")), 
       (gr_float *) field_array(field_id("energy")),
       (gr_float *) field_array(field_id("HI_density")), 
       (gr_float *) field_array(field_id("HII_density")), 
       (gr_float *) field_array(field_id("HM_density")),
       (gr_float *) field_array(field_id("HeI_density")), 
       (gr_float *) field_array(field_id("HeII_density")), 
       (gr_float *) field_array(field_id("HeIII_density")),
       (gr_float *) field_array(field_id("H2I_density")), 
       (gr_float *) field_array(field_id("H2II_density")),
       (gr_float *) field_array(field_id("DI_density")), 
       (gr_float *) field_array(field_id("DII_density")), 
       (gr_float *) field_array(field_id("HDI_density")),
       (gr_float *) field_array(field_id("e_density")), 
       (gr_float *) field_array(field_id("metal_density")),
       (gr_float *) field_array(field_id("pressure"))) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in calculate_pressure.\n");
  }

  if (calculate_gamma
      (chemistry_, units_,
       rank, m,
       (gr_float *) field_array(field_id("density")),
       (gr_float *) field_array(field_id("energy")),
       (gr_float *) field_array(field_id("HI_density")),
       (gr_float *) field_array(field_id("HII_density")),
       (gr_float *) field_array(field_id("HM_density")),
       (gr_float *) field_array(field_id("HeI_density")),
       (gr_float *) field_array(field_id("HeII_density")),
       (gr_float *) field_array(field_id("HeIII_density")),
       (gr_float *) field_array(field_id("H2I_density")),
       (gr_float *) field_array(field_id("H2II_density")),
       (gr_float *) field_array(field_id("DI_density")),
       (gr_float *) field_array(field_id("DII_density")),
       (gr_float *) field_array(field_id("HDI_density")),
       (gr_float *) field_array(field_id("e_density")),
       (gr_float *) field_array(field_id("metal_density")),
       (gr_float *) field_array(field_id("gamma"))) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in calculate_gamma.\n");
  }

#endif /* CONFIG_USE_GRACKLE */

}

//----------------------------------------------------------------------

double EnzoMethodGrackle::timestep ( CommBlock * comm_block ) throw()
{
#ifdef CONFIG_USE_GRACKLE
  initialize_(comm_block);
  return std::numeric_limits<double>::max();
#else
  return 0.0;
#endif /* CONFIG_USE_GRACKLE */
}

//======================================================================

