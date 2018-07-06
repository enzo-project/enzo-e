// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodGrackle class

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------------

EnzoMethodGrackle::EnzoMethodGrackle 
(
  EnzoConfig * config,
  const FieldDescr * field_descr
)
  : Method()
//#ifdef CONFIG_USE_GRACKLE
//  , chemistry_(0),
//    units_(0)
//#endif /* CONFIG_USE_GRACKLE */

{
#ifdef CONFIG_USE_GRACKLE

  /// Initialize default Refresh
  int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
		       enzo_sync_id_method_grackle);
  refresh(ir)->add_all_fields(field_descr->field_count());

  /// Initialize parameters

  units_     = config->method_grackle_units; // references?
  chemistry_ = config->method_grackle_chemistry;

  const gr_float a_value = 
    1. / (1. + config->physics_cosmology_initial_redshift);

  
  printf ("TRACE %s:%d calling initialize_chemistry_data\n",__FILE__,__LINE__);

  if (initialize_chemistry_data(&units_) == 0) {
      ERROR("EnzoMethodGrackle::EnzoMethodGrackle()",
      "Error in initialize_chemistry_data");
  }

  if (set_default_chemistry_parameters(chemistry_) == 0) {
      ERROR("EnzoMethodGrackle::EnzoMethodGrackle()",
      "Error in set_default_chemistry_parameters");
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

  // p | *chemistry_;
  WARNING ("EnzoMethodGrackle::pup()",
     "p | *chemistry_ not called!");
  //  p | *units_;
  WARNING ("EnzoMethodGrackle::pup()",
     "p | *units_ not called!");

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

  //  initialize_(block);

  Field field = block->data()->field();

  // Setup Grackle field struct for storing field data
  grackle_field_data grackle_fields_;

  // ASSUMES ALL ARRAYS ARE THE SAME SIZE

  gr_int m[3];
  field.dimensions (0,(int*)m,(int*)m+1,(int*)m+2);

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.field_size (0,&nx,&ny,&nz);

  gr_int grid_start[3] = {gx,      gy,      gz};
  gr_int grid_end[3]   = {gx+nx-1, gy+ny-1, gz+nz-1} ;

  // ASSUMES COSMOLOGY = false
  double a_value = 1.0;

  gr_int rank = block->rank();

  // Grackle grid dimenstion and grid size

  grackle_fields_.grid_rank      = rank;
  grackle_fields_.grid_dimension = new int[3];
  grackle_fields_.grid_start     = new int[3];
  grackle_fields_.grid_end       = new int[3];

  for (int i=0; i<3; i++){
    grackle_fields_.grid_dimension[i] = 1;
    grackle_fields_.grid_start[i]     = grid_start[i];
    grackle_fields_.grid_end[i]       = grid_end[i];
  }

  //grackle_fields_.grid_dimension[0] = nx;
  //grackle_fields_.grid_end[0] = nx - 1;

  // Setup all fields

  grackle_fields_.density         = (gr_float *) field.values("density");
  grackle_fields_.internal_energy = (gr_float *) field.values("internal_energy");
  grackle_fields_.x_velocity      = (gr_float *) field.values("velocity_x");
  grackle_fields_.y_velocity      = (gr_float *) field.values("velocity_y");
  grackle_fields_.z_velocity      = (gr_float *) field.values("velocity_z");
  grackle_fields_.HI_density      = (gr_float *) field.values("HI_density");
  grackle_fields_.HII_density     = (gr_float *) field.values("HII_density");
  grackle_fields_.HM_density      = (gr_float *) field.values("HM_density");
  grackle_fields_.HeI_density     = (gr_float *) field.values("HeI_density");
  grackle_fields_.HeII_density    = (gr_float *) field.values("HeII_density");
  grackle_fields_.HeIII_density   = (gr_float *) field.values("HeIII_density");
  grackle_fields_.H2I_density     = (gr_float *) field.values("H2I_density");
  grackle_fields_.H2II_density    = (gr_float *) field.values("H2II_density");
  grackle_fields_.DI_density      = (gr_float *) field.values("DI_density");
  grackle_fields_.DII_density     = (gr_float *) field.values("DII_density");
  grackle_fields_.HDI_density     = (gr_float *) field.values("HDI_density");
  grackle_fields_.e_density       = (gr_float *) field.values("e_density");
  grackle_fields_.metal_density   = (gr_float *) field.values("metal_density");
  //grackle_fields_.cooling_time  = (gr_float *) field.values("cooling_time");
  //grackle_fields_.temperature   = (gr_float *) field.values("temperature");
  //grackle_fields_.pressure      = (gr_float *) field.values("pressure");
  //grackle_fields_.gamma         = (gr_float *) field.values("gamma");

  double dt = block->dt();
  // double dt = 3.15e7 * 1e6 / my_units.time_units;

  if (solve_chemistry(&units_, &grackle_fields_, dt) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in solve_chemistry.\n");
  }

  // maybe change to "gr_float * cooling_time = grackle_fields_.cooling_time"
  gr_float *cooling_time;
  cooling_time = new gr_float[nx];

  if (calculate_cooling_time(&units_, &grackle_fields_, cooling_time) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in calculate_cooling_time.\n");
  }

  gr_float *temperature;
  temperature = new gr_float[nx];

  if (calculate_temperature(&units_, &grackle_fields_, temperature) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in calculate_temperature.\n");
  }

  gr_float *pressure;
  pressure = new gr_float[nx];

  if (calculate_pressure(&units_, &grackle_fields_, pressure) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in calculate_pressure.\n");
  }

  gr_float *gamma;
  gamma = new gr_float[nx];

  if (calculate_gamma(&units_, &grackle_fields_, gamma) == 0) {
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
  return std::numeric_limits<double>::max();
#else
  return 0.0;
#endif /* CONFIG_USE_GRACKLE */
}

//======================================================================


