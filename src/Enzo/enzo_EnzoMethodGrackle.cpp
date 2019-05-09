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

  // Gather list of fields that MUST be defined for this
  // method and check that they are permanent. If not,
  // define them:

  std::vector<std::string> fields_to_define;

  if (grackle_data->metal_cooling){

    if (! (field_descr->is_field("metal_density"))){
      fields_to_define.push_back("metal_density");

//      WARNING("EnzoMethodGrackle: ",
//              "Must define metal_density to use metal cooling with Grackle. Defining");
    }
  }

  if (grackle_data->primordial_chemistry > 0){
    std::string pc1_fields[6] = {"HI_density","HII_density",
                                  "HeI_density","HeII_density","HeIII_density",
                                  "e_density"};
    int numfields = 6;

    for(int ifield = 0; ifield < numfields; ifield++){
      if (! (field_descr->is_field( pc1_fields[ifield] ))){
        fields_to_define.push_back( pc1_fields[ifield] );

//        WARNING("EnzoMethodGrackle: ",
//                "Must define " + pc1_fields[ifield] + "if using primordial_chemistry = 1 with Grackle. Defining");
      }
    }

    if(grackle_data->primordial_chemistry > 1){

      std::string pc2_fields[3] = {"HM_density", "H2I_density", "H2II_density"};
      numfields = 3;

      for (int ifield = 0; ifield < numfields; ifield++){
        if (! (field_descr->is_field( pc2_fields[ifield] ))){
          fields_to_define.push_back( pc2_fields[ifield] );

//          WARNING("EnzoMethodGrackle: ",
//                  "Must define " + pc2_fields[ifield] + "if using primordial_chemistry = 2 with Grackle. Defining");
        }
      }

      if(grackle_data->primordial_chemistry > 2){
        std::string pc3_fields[3] = {"DI_density", "DII_density", "HDI_density"};
        numfields = 3;

        for(int ifield = 0; ifield < numfields; ifield++){
          if (! (field_descr->is_field( pc3_fields[ifield] ))){
            fields_to_define.push_back( pc3_fields[ifield] );

//            WARNING("EnzoMethodGrackle: ",
 //                   "Must define " + pc3_fields[ifield] +  "if using primordial_chemistry = 3 with Grackle. Defining");
          }
        }

      } // endif primordial_chemistry > 2

    } // endif primordial_chemistry > 1

  } // endif primordial chemistry is on

  if (grackle_data->use_specific_heating_rate){

    if ( !(field_descr->is_field("specific_heating_rate"))){
      fields_to_define.push_back("specific_heating_rate");

//      WARNING("EnzoMethodGrackle : ",
//             " Must define specific_heating_rate if using specific heating rate with Grackle. Defining")
    }
  }

  if (grackle_data->use_volumetric_heating_rate){
    if ( !(field_descr->is_field("volumetric_heating_rate"))){
      fields_to_define.push_back("volumetric_heating_rate");

//      WARNING("EnzoMethodGrackle: ",
//             "Must define volumetric_heating_rate field if using volumetric heating with Grackle. Defining.");
    }
  }

  // now define the fields
  for (int ifield = 0; ifield < fields_to_define.size(); ifield++){
    field_descr->insert_permanent( fields_to_define[ifield] );
  }

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

    EnzoPhysicsCosmology * cosmology = enzo::cosmology();

    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dt,
                                        time);
    grackle_units_.a_units
         = 1.0 / (1.0 + physics_cosmology_initial_redshift);
    grackle_units_.a_value = cosmo_a;

  } else if (enzo_config->method_grackle_radiation_redshift > -1){
    grackle_units_.a_value = 1.0 / (1.0 + enzo_config->method_grackle_radiation_redshift);
  }

  // Initialize grackle units and data
  TRACE("Calling initialize_chemistry_data from EnzoMethodGrackle::EnzoMethodGrackle()");
  if (initialize_chemistry_data(&grackle_units_) == ENZO_FAIL) {
    ERROR("EnzoConfig::EnzoConfig()",
    "Error in initialize_chemistry_data");
  }

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
    EnzoBlock * enzo_block = enzo::block(block);

    this->compute_(enzo_block);

    enzo_block->compute_done();
    return;
  #endif

}

#ifdef CONFIG_USE_GRACKLE
void EnzoMethodGrackle::setup_grackle_units (EnzoBlock * enzo_block,
                                             code_units * grackle_units,
                                             int i_hist /* default 0 */
                                             ) throw()
{
  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();

  /* Set code units for use in grackle */

  grackle_units->comoving_coordinates = enzo_config->physics_cosmology;
  grackle_units->density_units = enzo_units->density();
  grackle_units->length_units  = enzo_units->length();
  grackle_units->time_units    = enzo_units->time();
  grackle_units->velocity_units = enzo_units->velocity();

  grackle_units->a_units       = 1.0;
  grackle_units->a_value       = 1.0;
  if (grackle_units->comoving_coordinates){
    enzo_float cosmo_a  = 1.0;
    enzo_float cosmo_dt = 0.0;

    EnzoPhysicsCosmology * cosmology = enzo::cosmology();

    enzo_float compute_time;
    if (i_hist == 0) {
      compute_time = enzo_block->time();
    } else {
      Field field = enzo_block->data()->field();
      compute_time = field.history_time(i_hist);
    }

    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dt,
                                        compute_time);
    grackle_units->a_units
         = 1.0 / (1.0 + enzo_config->physics_cosmology_initial_redshift);
    grackle_units->a_value = cosmo_a;

  } else if (enzo_config->method_grackle_radiation_redshift > -1){
    grackle_units->a_value = 1.0 / (1.0 + enzo_config->method_grackle_radiation_redshift);
  }

  return;
}

void EnzoMethodGrackle::setup_grackle_fields(EnzoBlock * enzo_block,
  // Setup Grackle field struct for storing field data that will be passed
  // into Grackle. Initialize fields, if true, will also assign values to
  // the fields (equal to uniform background values). This is meant to be
  // used at initialization, and is by default false.
                                             grackle_field_data * grackle_fields_,
                                             int i_hist /*default 0 */
                                             ) throw()
  {

  Field field = enzo_block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  gr_int grid_dimension[3] = {ngx, ngy, ngz};
  gr_int grid_start[3]     = {gx,      gy,      gz};
  gr_int grid_end[3]       = {gx+nx-1, gy+ny-1, gz+nz-1} ;

  const gr_int rank = cello::rank();

  // Grackle grid dimenstion and grid size
  grackle_fields_->grid_rank      = rank;
  grackle_fields_->grid_dimension = new int[3];
  grackle_fields_->grid_start     = new int[3];
  grackle_fields_->grid_end       = new int[3];

  for (int i=0; i<3; i++){
    grackle_fields_->grid_dimension[i] = grid_dimension[i];
    grackle_fields_->grid_start[i]     = grid_start[i];
    grackle_fields_->grid_end[i]       = grid_end[i];
  }

  double hx, hy, hz;
  enzo_block->cell_width(&hx,&hy,&hz);
  grackle_fields_->grid_dx = hx;

  // Setup all fields to be passed into grackle
  grackle_fields_->density         = (gr_float *) field.values("density", i_hist);
  grackle_fields_->internal_energy = (gr_float *) field.values("internal_energy", i_hist);
  grackle_fields_->x_velocity      = (gr_float *) field.values("velocity_x", i_hist);
  grackle_fields_->y_velocity      = (gr_float *) field.values("velocity_y", i_hist);
  grackle_fields_->z_velocity      = (gr_float *) field.values("velocity_z", i_hist);

  // Get chemical species fields if they exist

  // primordial_chemistry == 0 fields
  grackle_fields_->HI_density      = field.is_field("HI_density") ?
                       (gr_float *) field.values("HI_density", i_hist)     : NULL;
  grackle_fields_->HII_density     = field.is_field("HII_density") ?
                       (gr_float *) field.values("HII_density", i_hist)    : NULL;;
  grackle_fields_->HeI_density     = field.is_field("HeI_density") ?
                       (gr_float *) field.values("HeI_density", i_hist)    : NULL;
  grackle_fields_->HeII_density    = field.is_field("HeII_density") ?
                       (gr_float *) field.values("HeII_density", i_hist)   : NULL;
  grackle_fields_->HeIII_density   = field.is_field("HeIII_density") ?
                       (gr_float *) field.values("HeIII_density", i_hist)  : NULL;
  grackle_fields_->e_density       = field.is_field("e_density") ?
                       (gr_float *) field.values("e_density", i_hist)      : NULL;

  // primordial_chemistry == 1 fields
  grackle_fields_->HM_density      = field.is_field("HM_density") ?
                       (gr_float *) field.values("HM_density", i_hist)     : NULL;
  grackle_fields_->H2I_density     = field.is_field("H2I_density") ?
                       (gr_float *) field.values("H2I_density", i_hist) : NULL;
  grackle_fields_->H2II_density    = field.is_field("H2II_density") ?
                       (gr_float *) field.values("H2II_density", i_hist) : NULL;

  // primordial_chemistry == 2 fields
  grackle_fields_->DI_density      = field.is_field("DI_density") ?
                       (gr_float *) field.values("DI_density", i_hist) : NULL;
  grackle_fields_->DII_density     = field.is_field("DII_density") ?
                       (gr_float *) field.values("DII_density", i_hist) : NULL;
  grackle_fields_->HDI_density     = field.is_field("HDI_density") ?
                       (gr_float *) field.values("HDI_density", i_hist) : NULL;

  grackle_fields_->metal_density   = field.is_field("metal_density") ?
                       (gr_float *) field.values("metal_density", i_hist) : NULL;

  /* Leave these as NULL for now and save for future development */
  gr_float * volumetric_heating_rate = NULL;
  gr_float * specific_heating_rate = NULL;

  grackle_fields_->volumetric_heating_rate = volumetric_heating_rate;
  grackle_fields_->specific_heating_rate   = specific_heating_rate;

  return;
}

void EnzoMethodGrackle::update_grackle_density_fields(
                               EnzoBlock * enzo_block,
                               grackle_field_data * grackle_fields_
                               ) throw() {

  // Intended for use at problem initialization. Scale species
  // density fields to be sensible mass fractions of the initial
  // density field. Problem types that require finer-tuned control
  // over individual species fields should adapt this function
  // in their initialization routines.

  Field field = enzo_block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  double tiny_number = 1.0E-10;

  for (int iz = 0; iz<ngz; iz++){
    for (int iy=0; iy<ngy; iy++){
      for (int ix=0; ix<ngx; ix++){
        int i = INDEX(ix,iy,iz,ngx,ngy);

        if(grackle_data->primordial_chemistry > 0){
          grackle_fields_->HI_density[i]   = grackle_fields_->density[i] * grackle_data->HydrogenFractionByMass;
          grackle_fields_->HII_density[i]   = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->HeI_density[i]   = grackle_fields_->density[i] * (1.0 - grackle_data->HydrogenFractionByMass);
          grackle_fields_->HeII_density[i]  = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->HeIII_density[i] = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->e_density[i]     = grackle_fields_->density[i] * tiny_number;
        }

        if (grackle_data->primordial_chemistry > 1){
          grackle_fields_->HM_density[i]    = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->H2I_density[i]   = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->H2II_density[i]  = grackle_fields_->density[i] * tiny_number;
        }

        if (grackle_data->primordial_chemistry > 2){
          grackle_fields_->DI_density[i]    = grackle_fields_->density[i] * grackle_data->DeuteriumToHydrogenRatio;
          grackle_fields_->DII_density[i]   = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->HDI_density[i]   = grackle_fields_->density[i] * tiny_number;
        }

      }
    }
  }

  return;
}

//----------------------------------------------------------------------
void EnzoMethodGrackle::compute_ ( EnzoBlock * enzo_block) throw()
{
  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();

  Field field = enzo_block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  const int rank = cello::rank();

  /* Set code units for use in grackle */
  grackle_field_data grackle_fields_;

  setup_grackle_units(enzo_block, &this->grackle_units_);
  setup_grackle_fields(enzo_block, &grackle_fields_);

  // Solve chemistry
  double dt = enzo_block->dt;
  if (solve_chemistry(&grackle_units_, &grackle_fields_, dt) == ENZO_FAIL) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in solve_chemistry.\n");
  }

  /* Correct total energy for changes in internal energy */
  enzo_float * total_energy    = (enzo_float *) field.values("total_energy");
  for (int i = 0; i < ngx*ngy*ngz; i++){
    total_energy[i] = grackle_fields_.internal_energy[i] +
            0.5 * grackle_fields_.x_velocity[i] * grackle_fields_.x_velocity[i];
    if (rank > 1) total_energy[i] += 0.5 * grackle_fields_.y_velocity[i] * grackle_fields_.y_velocity[i];
    if (rank > 2) total_energy[i] += 0.5 * grackle_fields_.z_velocity[i] * grackle_fields_.z_velocity[i];
  }

  // For testing purposes - reset internal energies with changes in mu
  if (enzo_config->initial_grackle_test_reset_energies){
    this->ResetEnergies(enzo_block);
  }

  delete_grackle_fields(&grackle_fields_);

  return;
}
#endif // config use grackle

//----------------------------------------------------------------------

double EnzoMethodGrackle::timestep ( Block * block ) throw()
{
  const EnzoConfig * config = enzo::config();

  double dt = std::numeric_limits<double>::max();;

#ifdef CONFIG_USE_GRACKLE
  if (config->method_grackle_use_cooling_timestep){
    EnzoBlock * enzo_block = enzo::block(block);
    Field field = enzo_block->data()->field();

    enzo_float * cooling_time = field.is_field("cooling_time") ?
                        (enzo_float *) field.values("cooling_time") : NULL;

    // make it if it doesn't exist
    bool delete_cooling_time = false;
    int gx,gy,gz;
    field.ghost_depth (0,&gx,&gy,&gz);

    int nx,ny,nz;
    field.size (&nx,&ny,&nz);

    int ngx = nx + 2*gx;
    int ngy = ny + 2*gy;
    int ngz = nz + 2*gz;

    int size = ngx*ngy*ngz;

    if (!(cooling_time)){
      cooling_time = new enzo_float [size];
      delete_cooling_time = true;
    }

    grackle_field_data grackle_fields_;

    setup_grackle_units(enzo_block,  &grackle_units_);
    setup_grackle_fields(enzo_block, &grackle_fields_);

    if (calculate_cooling_time(&grackle_units_, &grackle_fields_, cooling_time) == ENZO_FAIL) {
      ERROR("EnzoMethodGrackle::compute()",
      "Error in calculate_cooling_time.\n");
    }

    for (int i = 0; i < size; i++) dt = std::min(dt, std::abs(cooling_time[i]));

    if (delete_cooling_time){
      delete [] cooling_time;
    }

    delete_grackle_fields(&grackle_fields_);


  }
#endif

  return dt;
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
                                (enzo_float*) field.values("HI_density")    : NULL;
   enzo_float * HII_density   = field.is_field("HII_density") ?
                                (enzo_float*) field.values("HII_density")   : NULL;
   enzo_float * HeI_density   = field.is_field("HeI_density") ?
                                (enzo_float*) field.values("HeI_density")   : NULL;
   enzo_float * HeII_density  = field.is_field("HeII_density") ?
                                (enzo_float*) field.values("HeII_density")  : NULL;
   enzo_float * HeIII_density = field.is_field("HeIII_density") ?
                                (enzo_float*) field.values("HeIII_density") : NULL;
   enzo_float * e_density     = field.is_field("e_density") ?
                                (enzo_float*) field.values("e_density")     : NULL;

   enzo_float * H2I_density   = field.is_field("H2I_density") ?
                                (enzo_float*) field.values("H2I_density")   : NULL;
   enzo_float * H2II_density  = field.is_field("H2II_density") ?
                                (enzo_float*) field.values("H2II_density")  : NULL;
   enzo_float * HM_density    = field.is_field("HM_density") ?
                                (enzo_float*) field.values("HM_density")    : NULL;

   enzo_float * DI_density    = field.is_field("DI_density") ?
                               (enzo_float *) field.values("DI_density")    : NULL;
   enzo_float * DII_density   = field.is_field("DII_density") ?
                               (enzo_float *) field.values("DII_density")   : NULL;
   enzo_float * HDI_density   = field.is_field("HDI_density") ?
                               (enzo_float *) field.values("HDI_density")   : NULL;


   enzo_float * metal_density = field.is_field("metal_density") ?
                                (enzo_float*) field.values("metal_density") : NULL;

   // Field size
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

   double temperature_slope = log10(enzo_config->initial_grackle_test_maximum_temperature/
                                    enzo_config->initial_grackle_test_minimum_temperature)/
                                    double(ny);

   for (int iz=gz; iz<nz+gz; iz++){ // Metallicity
     for (int iy=gy; iy<ny+gy; iy++) { // Temperature
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

         internal_energy[i] = pow(10.0, ((temperature_slope * (iy-gy)) +
                              log10(enzo_config->initial_grackle_test_minimum_temperature)))/
                              mu / enzo_units->temperature() / (enzo_config->field_gamma - 1.0);
         total_energy[i] = internal_energy[i];

       }
     }
   }

  return;
}
#endif CONFIG_USE_GRACKLE

//======================================================================
