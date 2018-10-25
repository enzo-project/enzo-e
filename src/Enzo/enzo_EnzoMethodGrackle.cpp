// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Fri Apr  2 17:05:23 PDT 2010
/// @brief    Implements the EnzoMethodGrackle class

#include "cello.hpp"
#include "enzo.hpp"
extern CProxy_EnzoSimulation proxy_enzo_simulation;

//----------------------------------------------------------------------------

EnzoMethodGrackle::EnzoMethodGrackle
(
//  const FieldDescr * field_descr,
  const float physics_cosmology_initial_redshift,
  const float time
)
  : Method()
{
#ifdef CONFIG_USE_GRACKLE

  FieldDescr * field_descr = cello::field_descr();
  
  /// Initialize default Refresh
  int ir = add_refresh(4,0,neighbor_leaf,sync_barrier,
		       enzo_sync_id_method_grackle);
  refresh(ir)->add_all_fields();

  if (grackle_data->metal_cooling){
    ASSERT("EnzoMethodGrackle:Must define metal_density field ",
           "to use metal cooling with Grackle",
           field_descr->is_field("metal_density"));
  }

  if (grackle_data->primordial_chemistry > 0){
    ASSERT("EnzoMethodGrackle:Must define ",
           " HI_density field if using primordial_chemistry = 1 with Grackle",
           field_descr->is_field("HI_density"));
    ASSERT("EnzoMethodGrackle:Must define ",
           " HII_density field if using primordial_chemistry = 1 with Grackle",
           field_descr->is_field("HII_density"));
    ASSERT("EnzoMethodGrackle:Must define ",
           " HeI_density field if using primordial_chemistry = 1 with Grackle",
           field_descr->is_field("HeI_density"));
    ASSERT("EnzoMethodGrackle:Must define ",
           " HeII_density field if using primordial_chemistry = 1 with Grackle",
           field_descr->is_field("HeII_density"));
    ASSERT("EnzoMethodGrackle:Must define ",
           " HeIII_density field if using primordial_chemistry = 1 with Grackle",
           field_descr->is_field("HeIII_density"));
    ASSERT("EnzoMethodGrackle:Must define ",
           " e_density field if using primordial_chemistry = 1 with Grackle",
           field_descr->is_field("e_density"));

    if(grackle_data->primordial_chemistry > 1){
      ASSERT("EnzoMethodGrackle:Must define ",
             " HM_density field if using primordial_chemistry = 2 with Grackle",
             field_descr->is_field("HM_density"));
      ASSERT("EnzoMethodGrackle:Must define ",
             " H2I_density field if using primordial_chemistry = 2 with Grackle",
             field_descr->is_field("H2I_density"));
      ASSERT("EnzoMethodGrackle:Must define ",
             " H2II_density field if using primordial_chemistry = 2 with Grackle",
             field_descr->is_field("H2II_density"));
      if(grackle_data->primordial_chemistry > 2){
        ASSERT("EnzoMethodGrackle:Must define ",
               " DI_density field if using primordial_chemistry = 3 with Grackle",
               field_descr->is_field("DI_density"));
        ASSERT("EnzoMethodGrackle:Must define ",
               " DII_density field if using primordial_chemistry = 3 with Grackle",
               field_descr->is_field("DII_density"));
        ASSERT("EnzoMethodGrackle:Must define ",
               " HDI_density field if using primordial_chemistry = 3 with Grackle",
               field_descr->is_field("HDI_density"));
      }

    }

  } // field checks if primordial chemistry is on

  if (grackle_data->use_specific_heating_rate){
    ASSERT("EnzoMethodGrackle:Must define specific_heating_rate",
           " if using specific heating rate with Grackle",
           field_descr->is_field("specific_heating_rate"));
  }

  if (grackle_data->use_volumetric_heating_rate){
    ASSERT("EnzoMethodGrackle:Must define volumetric_heating_rate field ",
           " if using volumetric heating with Grackle",
           field_descr->is_field("volumetric_heating_rate"));
  }

  EnzoSimulation * simulation = proxy_enzo_simulation.ckLocalBranch();
  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();

  grackle_units_.comoving_coordinates = enzo_config->physics_cosmology;
  // Copy over code units to grackle units struct
  grackle_units_.density_units = enzo_units->density();
  grackle_units_.length_units  = enzo_units->length();
  grackle_units_.time_units    = enzo_units->time();
  grackle_units_.velocity_units = enzo_units->velocity();


  grackle_units_.a_units       = 1.0;
  grackle_units_.a_value       = 1.0;

  if (grackle_units_.comoving_coordinates){
    enzo_float cosmo_a  = 1.0;
    enzo_float cosmo_dt = 0.0;

    EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology *)
                simulation->problem()->physics("cosmology");
    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dt,
                                        time);
    grackle_units_.a_units
         = 1.0 / (1.0 + physics_cosmology_initial_redshift);
    grackle_units_.a_value = cosmo_a;

  }

  if (initialize_chemistry_data(&grackle_units_) == 0) {
    ERROR("EnzoConfig::EnzoConfig()",
    "Error in initialize_chemistry_data");
  }
  // method_grackle_units = grackle_units_; // copy over to global
  printf ("TRACE %s:%d calling initialize_chemistry_data\n",__FILE__,__LINE__);

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

// void EnzoMethodGrackle::pup (PUP::er &p) - in hpp


void EnzoMethodGrackle::compute ( Block * block) throw()
{
  if (!block->is_leaf()) return;

  #ifndef CONFIG_USE_GRACKLE

    ERROR("EnzoMethodGrackle::compute()",
    "Trying to use method 'grackle' with "
    "Grackle configuration turned off!");

  #else /* CONFIG_USE_GRACKLE */
    EnzoBlock * enzo_block = static_cast<EnzoBlock*> (block);

    this->compute_(enzo_block);

    enzo_block->compute_done();
    return;
  #endif

}

//----------------------------------------------------------------------
#ifdef CONFIG_USE_GRACKLE
void EnzoMethodGrackle::compute_ ( EnzoBlock * enzo_block) throw()
{
  EnzoSimulation * simulation = proxy_enzo_simulation.ckLocalBranch();
  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();

  /* Set code units for use in grackle */

  grackle_units_.comoving_coordinates = enzo_config->physics_cosmology;
  grackle_units_.density_units = enzo_units->density();
  grackle_units_.length_units  = enzo_units->length();
  grackle_units_.time_units    = enzo_units->time();
  grackle_units_.velocity_units = enzo_units->velocity();

  grackle_units_.a_units       = 1.0;
  grackle_units_.a_value       = 1.0;
  if (grackle_units_.comoving_coordinates){
    enzo_float cosmo_a  = 1.0;
    enzo_float cosmo_dt = 0.0;

    EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology *)
                simulation->problem()->physics("cosmology");
    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dt,
                                        enzo_block->time());
    grackle_units_.a_units
         = 1.0 / (1.0 + enzo_config->physics_cosmology_initial_redshift);
    grackle_units_.a_value = cosmo_a;

  }

  //  initialize_(block);
  Field field = enzo_block->data()->field();

  // Setup Grackle field struct for storing field data
  grackle_field_data grackle_fields_;

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.field_size (0,&nx,&ny,&nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  gr_int grid_dimension[3] = {ngx, ngy, ngz};
  gr_int grid_start[3]     = {gx,      gy,      gz};
  gr_int grid_end[3]       = {gx+nx-1, gy+ny-1, gz+nz-1} ;

/*
  EnzoPhysicsCosmology * cosmology = (EnzoPhysicsCosmology *)
    simulation->problem()->physics("cosmology");

  if (grackle_units_.comoving_coordinates){
    enzo_float cosmo_a = 1.0, cosmo_dadt = 0.0;

    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dadt,
                                        enzo_block->time());
    grackle_units_.a_value = cosmo_a;
  }
*/
  const gr_int rank = cello::rank();

  // Grackle grid dimenstion and grid size

  grackle_fields_.grid_rank      = rank;
  grackle_fields_.grid_dimension = new int[3];
  grackle_fields_.grid_start     = new int[3];
  grackle_fields_.grid_end       = new int[3];

  for (int i=0; i<3; i++){
    grackle_fields_.grid_dimension[i] = grid_dimension[i];
    grackle_fields_.grid_start[i]     = grid_start[i];
    grackle_fields_.grid_end[i]       = grid_end[i];
  }

  double hx, hy, hz;
  enzo_block->cell_width(&hx,&hy,&hz);
  grackle_fields_.grid_dx = hx;

  // Setup all fields to be passed into grackle
  grackle_fields_.density         = (gr_float *) field.values("density");
  grackle_fields_.internal_energy = (gr_float *) field.values("internal_energy");
  grackle_fields_.x_velocity      = (gr_float *) field.values("velocity_x");
  grackle_fields_.y_velocity      = (gr_float *) field.values("velocity_y");
  grackle_fields_.z_velocity      = (gr_float *) field.values("velocity_z");

  grackle_fields_.HI_density      = field.is_field("HI_density") ?
                       (gr_float *) field.values("HI_density")     : NULL;
  grackle_fields_.HII_density     = field.is_field("HII_density") ?
                       (gr_float *) field.values("HII_density")    : NULL;
  grackle_fields_.HM_density      = field.is_field("HM_density") ?
                       (gr_float *) field.values("HM_density")     : NULL;
  grackle_fields_.HeI_density     = field.is_field("HeI_density") ?
                       (gr_float *) field.values("HeI_density")    : NULL;
  grackle_fields_.HeII_density    = field.is_field("HeII_density") ?
                       (gr_float *) field.values("HeII_density")   : NULL;
  grackle_fields_.HeIII_density   = field.is_field("HeIII_density") ?
                       (gr_float *) field.values("HeIII_density")  : NULL;
  grackle_fields_.e_density       = field.is_field("e_density") ?
                       (gr_float *) field.values("e_density")      : NULL;


  grackle_fields_.H2I_density     = field.is_field("H2I_density") ?
                       (gr_float *) field.values("H2I_density") : NULL;
  grackle_fields_.H2II_density    = field.is_field("H2II_density") ?
                       (gr_float *) field.values("H2II_density") : NULL;
  grackle_fields_.DI_density      = field.is_field("DI_density") ?
                       (gr_float *) field.values("DI_density") : NULL;
  grackle_fields_.DII_density     = field.is_field("DII_density") ?
                       (gr_float *) field.values("DII_density") : NULL;
  grackle_fields_.HDI_density     = field.is_field("HDI_Density") ?
                       (gr_float *) field.values("HDI_density") : NULL;
  grackle_fields_.metal_density   = field.is_field("metal_density") ?
                       (gr_float *) field.values("metal_density") : NULL;

  gr_float * volumetric_heating_rate = NULL;
  gr_float * specific_heating_rate = NULL;

  grackle_fields_.volumetric_heating_rate = volumetric_heating_rate;
  grackle_fields_.specific_heating_rate   = specific_heating_rate;

  double dt = enzo_block->dt;

  if (solve_chemistry(&grackle_units_, &grackle_fields_, dt) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in solve_chemistry.\n");
  }

  /* if this is a test problem, reset energies */
  // AE: Need to fix this to only call when in test
  if (enzo_config->initial_grackle_test_reset_energies){
    this->ResetEnergies(enzo_block);
  }

  /* might want to not do it this way */
  const int in = cello::index_static();
  int comoving_coordinates = enzo_config->physics_cosmology;
  enzo_float * pressure    = field.is_field("pressure") ?
               (enzo_float*) field.values("pressure") : NULL;
  enzo_float * temperature = field.is_field("temperature") ?
               (enzo_float*) field.values("temperature") : NULL;
  if (pressure){
    EnzoComputePressure compute_pressure (EnzoBlock::Gamma[in],
                                          comoving_coordinates);
    compute_pressure.compute(enzo_block);
  }

  if (temperature){
    EnzoComputeTemperature compute_temperature
      (enzo_config->ppm_density_floor,
       enzo_config->ppm_temperature_floor,
       enzo_config->ppm_mol_weight,
       comoving_coordinates);

    compute_temperature.compute(enzo_block);
  }

  gr_float * cooling_time = field.is_field("cooling_time") ?
                    (gr_float *) field.values("cooling_time") : NULL;
  if (cooling_time){
    if (calculate_cooling_time(&grackle_units_, &grackle_fields_, cooling_time) == 0) {
      ERROR("EnzoMethodGrackle::compute()",
      "Error in calculate_cooling_time.\n");
    }
 }
/*
  gr_float *temperature;
  temperature = new gr_float[nx];

  if (calculate_temperature(&grackle_units_, &grackle_fields_, temperature) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in calculate_temperature.\n");
  }

  gr_float *pressure;
  pressure = new gr_float[nx];

  if (calculate_pressure(&grackle_units_, &grackle_fields_, pressure) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in calculate_pressure.\n");
  }

  gr_float *gamma;
  gamma = new gr_float[nx];

  if (calculate_gamma(&units_, &grackle_fields_, gamma) == 0) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in calculate_gamma.\n");
  }
  */

  return;
}
#endif // config use grackle

//----------------------------------------------------------------------

double EnzoMethodGrackle::timestep ( Block * block ) const throw()
{
#ifdef CONFIG_USE_GRACKLE
  // return std::numeric_limits<double>::max();
  return 1.0E-6;
#else
  return 0.0;
#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE
void EnzoMethodGrackle::ResetEnergies ( EnzoBlock * enzo_block) throw()
{
   const EnzoConfig * enzo_config = enzo::config();
   EnzoUnits * enzo_units = enzo::units();

   /* Only need to do this if tracking chemistry */
   if (grackle_data->primordial_chemistry < 1)
     return;

   Field field = enzo_block->data()->field();

   enzo_float * density     = (enzo_float*) field.values("density");
   enzo_float * internal_energy = (enzo_float*) field.values("internal_energy");
   enzo_float * total_energy    = (enzo_float*) field.values("total_energy");

   enzo_float * pressure    = field.is_field("pressure") ?
                (enzo_float*) field.values("pressure") : NULL;
   enzo_float * temperature = field.is_field("temperature") ?
                (enzo_float*) field.values("temperature") : NULL;

   enzo_float * HI_density    = field.is_field("HI_density") ?
                  (enzo_float*) field.values("HI_density") : NULL;
   enzo_float * HII_density   = field.is_field("HII_density") ?
                  (enzo_float*) field.values("HII_density") : NULL;
   enzo_float * HeI_density   = field.is_field("HeI_density") ?
                  (enzo_float*) field.values("HeI_density") : NULL;
   enzo_float * HeII_density  = field.is_field("H2II_density") ?
                  (enzo_float*) field.values("HeII_density") : NULL;
   enzo_float * HeIII_density = field.is_field("HeIII_density") ?
                                (enzo_float*) field.values("HeIII_density") : NULL;
   enzo_float * H2I_density   = field.is_field("H2I_density") ?
                                  (enzo_float*) field.values("H2I_density") : NULL;
   enzo_float * H2II_density  = field.is_field("H2II_density") ?
                                  (enzo_float*) field.values("H2II_density") : NULL;
   enzo_float * HM_density    = field.is_field("HM_density") ?
                                  (enzo_float*) field.values("HM_density") : NULL;
   enzo_float * DI_density      = field.is_field("DI_density") ?
                                  (enzo_float *) field.values("DI_density") : NULL;
   enzo_float * DII_density     = field.is_field("DII_density") ?
                                  (enzo_float *) field.values("DII_density") : NULL;
   enzo_float * HDI_density     = field.is_field("HDI_density") ?
                                  (enzo_float *) field.values("HDI_density"): NULL;
   enzo_float * e_density     = field.is_field("e_density") ?
                  (enzo_float*) field.values("e_density") : NULL;
   enzo_float * metal_density = field.is_field("metal_density") ?
                            (enzo_float*) field.values("metal_density") : NULL;

   int nx,ny,nz;
   field.size(&nx,&ny,&nz);

   // Cell widths
   double xm,ym,zm;
   enzo_block->data()->lower(&xm,&ym,&zm);
   double xp,yp,zp;
   enzo_block->data()->upper(&xp,&yp,&zp);

   // Ghost depths
   int gx,gy,gz;
   field.ghost_depth(0,&gx,&gy,&gz);

   int mx,my,mz;
   field.dimensions(0,&mx,&my,&mz);

   const int m = mx*my*mz;

   int ngx = nx + 2*gx;
   int ngy = ny + 2*gy;
   int ngz = nz + 2*gz;


   double a_units = 1.0 / (1.0 + enzo_config->physics_cosmology_initial_redshift);

   const double mh = 1.67262171E-24;
   const double kboltz = 1.3806504E-16;

   gr_float temperature_units =  mh * pow(a_units *
                                          enzo_units->velocity(), 2) / kboltz;

   double temperature_slope = log10(enzo_config->initial_grackle_test_maximum_temperature/
                                    enzo_config->initial_grackle_test_minimum_temperature)/
                                    double(nz);

   for (int iz=gz; iz<nz+gz; iz++){ // Temperature
     for (int iy=gy; iy<ny+gy; iy++) { // Metallicity
       for (int ix=gx; ix<nx+gx; ix++) { // H Number Density
         int i = INDEX(ix,iy,iz,ngx,ngy);

         enzo_float mu = e_density[i] + HI_density[i] + HII_density[i] +
            (HeI_density[i] + HeII_density[i] + HeIII_density[i])*0.25;

         if (grackle_data->primordial_chemistry > 1){
           mu += HM_density[i] + 0.5 * (H2I_density[i] + H2II_density[i]);
         }

         if (grackle_data->primordial_chemistry > 2){
           mu += (DI_density[i] + DII_density[i])*0.5 + HDI_density[i]/3.0;
         }

         mu = density[i] / mu;

         internal_energy[i] = pow(10.0, ((temperature_slope * (iz-gz)) +
                              log10(enzo_config->initial_grackle_test_minimum_temperature)))/
                              mu / temperature_units / (enzo_config->field_gamma - 1.0);
         total_energy[i] = internal_energy[i];

       }
     }
   }

  return;
}
#endif CONFIG_USE_GRACKLE

//======================================================================
