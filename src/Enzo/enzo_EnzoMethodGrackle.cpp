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
  , chemistry_(0),
    units_(0)
#endif /* CONFIG_USE_GRACKLE */

{
#ifdef CONFIG_USE_GRACKLE

  set_num_refresh(1);

  refresh(0)->set_ghost_depth(4);
  refresh(0)->set_min_face_rank(0);
  refresh(0)->add_all_fields(field_descr->field_count());
  refresh(0)->set_sync_type(sync_neighbor);

  /// Initialize parameters

  units_     = & config->method_grackle_units;
  chemistry_ = & config->method_grackle_chemistry;

  const gr_float a_value = 
    1. / (1. + config->physics_cosmology_initial_redshift);

  
  printf ("TRACE %s:%d calling initialize_chemistry_data\n");

  if (initialize_chemistry_data
      (*chemistry_, *units_, a_value) == 0) {
    ERROR("EnzoMethodGrackle::EnzoMethodGrackle()",
	  "Error in initialize_chemistry_data");
  }
  
#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

void EnzoMethodGrackle::pup (PUP::er &p)
{

  // NOTE: change this function whenever attributes change

#ifdef CONFIG_USE_GRACKLE

  TRACEPUP;

  Method::pup(p);

  p | *chemistry_;
  p | *units_;

#endif /* CONFIG_USE_GRACKLE */

}

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute ( Block * block) throw()
{

  if (!block->is_leaf()) return;

#ifndef CONFIG_USE_GRACKLE

  ERROR("EnzoMethodGrackle::compute()",
	"Trying to use method 'grackle' with "
	"Grackle configuration turned off!");

#else /* CONFIG_USE_GRACKLE */

  initialize_(block);

  Field field = block->data()->field();

  // ASSUMES ALL ARRAYS ARE THE SAME SIZE
  gr_int m[3];
  array_dimension (0,m,m+1,m+2);

  int gx,gy,gz;
  ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  block_size (&nx,&ny,&nz);

  gr_int grid_start[3] = {gx,      gy,      gz};
  gr_int grid_end[3]   = {gx+nx-1, gy+ny-1, gz+nz-1} ;

  // ASSUMES COSMOLOGY = false
  double a_value = 1.0;

  gr_int rank = this->rank();

  gr_float * density       = (gr_float *) field.values("density");
  gr_float * energy        = (gr_float *) field.values("internal_energy");
  gr_float * velocity_x    = (gr_float *) field.values("velocity_x");
  gr_float * velocity_y    = (gr_float *) field.values("velocity_y");
  gr_float * velocity_z    = (gr_float *) field.values("velocity_z");
  gr_float * HI_density    = (gr_float *) field.values("HI_density");
  gr_float * HII_density   = (gr_float *) field.values("HII_density");
  gr_float * HM_density    = (gr_float *) field.values("HM_density");
  gr_float * HeI_density   = (gr_float *) field.values("HeI_density");
  gr_float * HeII_density  = (gr_float *) field.values("HeII_density");
  gr_float * HeIII_density = (gr_float *) field.values("HeIII_density");
  gr_float * H2I_density   = (gr_float *) field.values("H2I_density");
  gr_float * H2II_density  = (gr_float *) field.values("H2II_density");
  gr_float * DI_density    = (gr_float *) field.values("DI_density");
  gr_float * DII_density   = (gr_float *) field.values("DII_density");
  gr_float * HDI_density   = (gr_float *) field.values("HDI_density");
  gr_float * e_density     = (gr_float *) field.values("e_density");
  gr_float * metal_density = (gr_float *) field.values("metal_density");
  gr_float * cooling_time  = (gr_float *) field.values("cooling_time");
  gr_float * temperature   = (gr_float *) field.values("temperature");
  gr_float * pressure      = (gr_float *) field.values("pressure");
  gr_float * gamma         = (gr_float *) field.values("gamma");

  double dt = time_step();

  if (solve_chemistry
      (*chemistry_, *units_,
       a_value, dt,
       rank, m,
       grid_start,  grid_end,
       density,     energy,
       velocity_x,  velocity_y,   velocity_z,
       HI_density,  HII_density,
       HM_density,
       HeI_density, HeII_density, HeIII_density,
       H2I_density, H2II_density,
       DI_density,  DII_density,
       HDI_density,
       e_density,
       metal_density) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in solve_chemistry");
  }

  if (calculate_cooling_time
      (*chemistry_, *units_,
       a_value,
       rank, m,
       grid_start,  grid_end,
       density,     energy,
       velocity_x,  velocity_y,   velocity_z,
       HI_density,  HII_density,
       HM_density,
       HeI_density, HeII_density, HeIII_density,
       H2I_density, H2II_density,
       DI_density,  DII_density,
       HDI_density,
       e_density,
       metal_density,
       cooling_time) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in calculate_cooling_time.\n");
  }

  if (calculate_temperature
      (*chemistry_, *units_,
       rank, m,
       density,     energy,
       HI_density,  HII_density,
       HM_density,
       HeI_density, HeII_density, HeIII_density,
       H2I_density, H2II_density,
       DI_density,  DII_density,
       HDI_density,
       e_density,
       metal_density,
       temperature) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in calculate_temperature.\n");
  }

  if (calculate_pressure
      (*chemistry_, *units_,
       rank, m,
       density,     energy,
       HI_density,  HII_density,
       HM_density,
       HeI_density, HeII_density,  HeIII_density,
       H2I_density, H2II_density,
       DI_density,  DII_density,
       HDI_density,
       e_density,
       metal_density,
       pressure) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in calculate_pressure.\n");
  }

  if (calculate_gamma
      (*chemistry_, *units_,
       rank, m,
       density,     energy,
       HI_density,  HII_density,
       HM_density,
       HeI_density, HeII_density,  HeIII_density,
       H2I_density, H2II_density,
       DI_density,  DII_density,
       HDI_density,
       e_density,
       metal_density,
       gamma) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
	  "Error in calculate_gamma.\n");
  }

#endif /* CONFIG_USE_GRACKLE */

  EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

  enzo_block->compute_done();

}

//----------------------------------------------------------------------

double EnzoMethodGrackle::timestep ( Block * block ) const throw()
{
#ifdef CONFIG_USE_GRACKLE
  initialize_(block);
  return std::numeric_limits<double>::max();
#else
  return 0.0;
#endif /* CONFIG_USE_GRACKLE */
}

//======================================================================

