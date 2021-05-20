// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
///           Andrew Emerick (aemerick11@gmail.com)
/// @date     Tues May  7
/// @brief    Implements the EnzoMethodGrackle class

#include "cello.hpp"
#include "enzo.hpp"


//----------------------------------------------------------------------------

EnzoMethodGrackle::EnzoMethodGrackle
(
  const double physics_cosmology_initial_redshift,
  const double time
)
  : Method()
#ifdef CONFIG_USE_GRACKLE
    , grackle_units_(),
    grackle_rates_(),
    time_grackle_data_initialized_(ENZO_FLOAT_UNDEFINED)
#endif
{
#ifdef CONFIG_USE_GRACKLE

  // Gather list of fields that MUST be defined for this
  // method and check that they are permanent. If not,
  // define them.
  // ----- EnzoMethodGrackle::define_required_grackle_fields();

  FieldDescr * field_descr = cello::field_descr();

  // special container for ensuring color fields are properly grouped
  const int rank = cello::rank();
  std::vector<std::string> color_fields;
  chemistry_data * grackle_chemistry =
      enzo::config()->method_grackle_chemistry;

  this->required_fields_ = std::vector<std::string> {"density","internal_energy",
                                                     "total_energy"};
  if (rank>=0) this->required_fields_.push_back("velocity_x");
  if (rank>=1) this->required_fields_.push_back("velocity_y");
  if (rank>=2) this->required_fields_.push_back("velocity_z");

  if (grackle_chemistry->metal_cooling > 0){
    this->required_fields_.push_back("metal_density");
    color_fields.push_back("metal_density");
  }

  // Define primordial chemistry fields
  if (grackle_chemistry->primordial_chemistry > 0){
    std::vector<std::string> pc1_fields {"HI_density","HII_density",
                                         "HeI_density","HeII_density","HeIII_density",
                                         "e_density"};

    this->required_fields_.insert(this->required_fields_.end(), pc1_fields.begin(), pc1_fields.end());
    color_fields.insert(color_fields.end(), pc1_fields.begin(), pc1_fields.end());

    if(grackle_chemistry->primordial_chemistry > 1){

      std::vector<std::string> pc2_fields {"HM_density", "H2I_density", "H2II_density"};
      this->required_fields_.insert(this->required_fields_.end(), pc2_fields.begin(), pc2_fields.end());
      color_fields.insert(color_fields.end(), pc2_fields.begin(), pc2_fields.end());

      if(grackle_chemistry->primordial_chemistry > 2){
        std::vector<std::string> pc3_fields {"DI_density", "DII_density", "HDI_density"};
        this->required_fields_.insert(this->required_fields_.end(), pc3_fields.begin(), pc3_fields.end());
        color_fields.insert(color_fields.end(), pc3_fields.begin(), pc3_fields.end());
      } // endif primordial_chemistry > 2
    } // endif primordial_chemistry > 1
  } // endif primordial chemistry is on

  if (grackle_chemistry->use_specific_heating_rate)
      this->required_fields_.push_back("specific_heating_rate");

  if (grackle_chemistry->use_volumetric_heating_rate)
      this->required_fields_.push_back("volumetric_heating_rate");

  // Define fields and assign fields to correct
  this->define_fields();
  this->define_group_fields(color_fields, "color");

  /// Initialize default Refresh
  cello::simulation()->new_refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  /// Define Grackle's internal data structures
  time_grackle_data_initialized_ = ENZO_FLOAT_UNDEFINED;
  initialize_grackle_chemistry_data(time);

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute ( Block * block) throw()
{

  if (block->is_leaf()){

  #ifndef CONFIG_USE_GRACKLE

    ERROR("EnzoMethodGrackle::compute()",
    "Trying to use method 'grackle' with "
    "Grackle configuration turned off!");

  #else /* CONFIG_USE_GRACKLE */

    EnzoBlock * enzo_block = enzo::block(block);

    // Start timer
    Simulation * simulation = cello::simulation();
    if (simulation)
      simulation->performance()->start_region(perf_grackle,__FILE__,__LINE__);

    this->compute_(enzo_block);

    enzo_block->compute_done();

    if (simulation)
      simulation->performance()->stop_region(perf_grackle,__FILE__,__LINE__);
  #endif
  }

  return;

}

#ifdef CONFIG_USE_GRACKLE

void EnzoMethodGrackle::define_required_grackle_fields()
{
  // Gather list of fields that MUST be defined for this method and
  // check that they are permanent. If not, define them.

  // This has been split off from the constructor so that other methods that
  // are initialized first and need knowledge of these fields at initialization
  // (e.g. to set up a refresh object), can ensure that the fields are defined

  if (!enzo::config()->method_grackle_use_grackle) {return;}

  FieldDescr * field_descr = cello::field_descr();
  Config   * config  = (Config *) cello::config();  

  // special container for ensuring color fields are properly grouped
  const int rank = cello::rank();
  std::vector<std::string> color_fields;
  chemistry_data * grackle_chemistry =
      enzo::config()->method_grackle_chemistry;

  std::vector<std::string> required_fields = std::vector<std::string> {"density","internal_energy",
                                                     "total_energy"};
  if (rank>=0) required_fields.push_back("velocity_x");
  if (rank>=1) required_fields.push_back("velocity_y");
  if (rank>=2) required_fields.push_back("velocity_z");

  if (grackle_chemistry->metal_cooling > 0){
    required_fields.push_back("metal_density");
    color_fields.push_back("metal_density");
  }

  // Define primordial chemistry fields
  if (grackle_chemistry->primordial_chemistry > 0){
    std::vector<std::string> pc1_fields {"HI_density","HII_density",
                                         "HeI_density","HeII_density","HeIII_density",
                                         "e_density"};

    required_fields.insert(required_fields.end(), pc1_fields.begin(), pc1_fields.end());
    color_fields.insert(color_fields.end(), pc1_fields.begin(), pc1_fields.end());

    if(grackle_chemistry->primordial_chemistry > 1){

      std::vector<std::string> pc2_fields {"HM_density", "H2I_density", "H2II_density"};
      required_fields.insert(required_fields.end(), pc2_fields.begin(), pc2_fields.end());
      color_fields.insert(color_fields.end(), pc2_fields.begin(), pc2_fields.end());

      if(grackle_chemistry->primordial_chemistry > 2){
        std::vector<std::string> pc3_fields {"DI_density", "DII_density", "HDI_density"};
        required_fields.insert(required_fields.end(), pc3_fields.begin(), pc3_fields.end());
        color_fields.insert(color_fields.end(), pc3_fields.begin(), pc3_fields.end());
      } // endif primordial_chemistry > 2
    } // endif primordial_chemistry > 1
  } // endif primordial chemistry is on

  if (grackle_chemistry->use_specific_heating_rate)
      required_fields.push_back("specific_heating_rate");

  if (grackle_chemistry->use_volumetric_heating_rate)
      required_fields.push_back("volumetric_heating_rate");

  // Define fields and assign fields to correct
  //this->define_fields();
  bool added_fields = false;

  for (int ifield = 0; ifield < required_fields.size(); ifield++){
    std::string field = required_fields[ifield];
    if( ! field_descr->is_field( field )){
      int id_field = field_descr->insert_permanent( field );

      field_descr->set_precision(id_field, config->field_precision);
      added_fields = true;
    }
  }

  //this->define_group_fields(color_fields, "color");
  for (int ifield = 0; ifield < color_fields.size(); ifield++){

    // Maybe just throw error here to keep this fully separate from above
    if( ! field_descr->is_field( required_fields[ifield] )){
      int field_id = field_descr->insert_permanent( required_fields[ifield] );
      field_descr->set_precision(field_id, config->field_precision);
      added_fields = true;
    }

    if (!(field_descr->groups()->is_in( color_fields[ifield], "color")) ){
      field_descr->groups()->add( color_fields[ifield], "color");
    }

  }

  // Need to reconstruct history if new fields added
  if (added_fields) field_descr->reset_history(config->field_history);


}

//----------------------------------------------------------------------

void EnzoMethodGrackle::initialize_grackle_chemistry_data(double current_time)
{

  /* Define Grackle's chemistry data if not yet defined */

  if (this->time_grackle_data_initialized_ == current_time) return;

  if (this->time_grackle_data_initialized_ != ENZO_FLOAT_UNDEFINED){
    // deallocate previously the allocated allocated grackle_rates_ (doesn't
    // actually affect the chemistry_data pointer)
    deallocate_grackle_rates_();
  }

  EnzoUnits * enzo_units = enzo::units();
  const EnzoConfig * enzo_config = enzo::config();

  grackle_units_.comoving_coordinates = enzo_config->physics_cosmology;

  // Copy over code units to grackle units struct
  grackle_units_.density_units  = enzo_units->density();
  grackle_units_.length_units   = enzo_units->length();
  grackle_units_.time_units     = enzo_units->time();
  grackle_units_.velocity_units = enzo_units->velocity();
  grackle_units_.a_units        = 1.0;
  grackle_units_.a_value        = 1.0;

  if (grackle_units_.comoving_coordinates){
    enzo_float cosmo_a  = 1.0;
    enzo_float cosmo_dt = 0.0;

    EnzoPhysicsCosmology * cosmology = enzo::cosmology();

    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dt,
                                        current_time);
    grackle_units_.a_units
         = 1.0 / (1.0 + enzo_config->physics_cosmology_initial_redshift);
    grackle_units_.a_value = cosmo_a;

  } else if (enzo_config->method_grackle_radiation_redshift > -1){
    grackle_units_.a_value = 1.0 /
                         (1.0 + enzo_config->method_grackle_radiation_redshift);
  }

  // Initialize grackle units and data
  TRACE("Calling initialize_chemistry_data from EnzoMethodGrackle::EnzoMethodGrackle()");

  if (_initialize_chemistry_data(enzo_config->method_grackle_chemistry,
				 &grackle_rates_, &grackle_units_)
      == ENZO_FAIL) {
    ERROR("EnzoMethodGrackle::initialize_grackle_chemistry_data",
    "Error in _initialize_chemistry_data");
  }

  this->time_grackle_data_initialized_ = current_time;

  return;
}

//----------------------------------------------------------------------------

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
    grackle_units->a_value = 1.0 /
                         (1.0 + enzo_config->method_grackle_radiation_redshift);
  }

  return;
}

//--------------------------------------------------------------------------

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

//----------------------------------------------------------------------------

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

  chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

  for (int iz = 0; iz<ngz; iz++){
    for (int iy=0; iy<ngy; iy++){
      for (int ix=0; ix<ngx; ix++){
        int i = INDEX(ix,iy,iz,ngx,ngy);

        if(grackle_chemistry->primordial_chemistry > 0){
          grackle_fields_->HI_density[i]   = grackle_fields_->density[i] * grackle_chemistry->HydrogenFractionByMass;
          grackle_fields_->HII_density[i]   = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->HeI_density[i]   = grackle_fields_->density[i] * (1.0 - grackle_chemistry->HydrogenFractionByMass);
          grackle_fields_->HeII_density[i]  = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->HeIII_density[i] = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->e_density[i]     = grackle_fields_->density[i] * tiny_number;
        }

        if (grackle_chemistry->primordial_chemistry > 1){
          grackle_fields_->HM_density[i]    = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->H2I_density[i]   = grackle_fields_->density[i] * tiny_number;
          grackle_fields_->H2II_density[i]  = grackle_fields_->density[i] * tiny_number;
        }

        if (grackle_chemistry->primordial_chemistry > 2){
          grackle_fields_->DI_density[i]    = grackle_fields_->density[i] * grackle_chemistry->DeuteriumToHydrogenRatio;
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

  chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

  // Solve chemistry
  double dt = enzo_block->dt;
  if (local_solve_chemistry(grackle_chemistry, &grackle_rates_,
			    &grackle_units_, &grackle_fields_, dt)
      == ENZO_FAIL) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in local_solve_chemistry.\n");
  }

  /* Correct total energy for changes in internal energy */
  gr_float * v3[3];
  v3[0] = grackle_fields_.x_velocity;
  v3[1] = grackle_fields_.y_velocity;
  v3[2] = grackle_fields_.z_velocity;

  const bool mhd = field.is_field("bfield_x");
  enzo_float * b3[3] = {NULL, NULL, NULL};
  if (mhd) {
    b3[0]                = (enzo_float*) field.values("bfield_x");
    if (rank >= 2) b3[1] = (enzo_float*) field.values("bfield_y");
    if (rank >= 3) b3[2] = (enzo_float*) field.values("bfield_z");
  }

  enzo_float * total_energy    = (enzo_float *) field.values("total_energy");
  for (int i = 0; i < ngx*ngy*ngz; i++){
    total_energy[i] = grackle_fields_.internal_energy[i];

    enzo_float inv_density;
    if (mhd) inv_density = 1.0 / grackle_fields_.density[i];
    for (int dim = 0; dim < rank; dim++){
      total_energy[i] += 0.5 * v3[dim][i] * v3[dim][i];
      if (mhd) total_energy[i] += 0.5 * b3[dim][i] * b3[dim][i] * inv_density;
    }
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

double EnzoMethodGrackle::timestep ( Block * block ) const throw()
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

    calculate_cooling_time(block, cooling_time, NULL, NULL, 0);

    // make sure to exclude the ghost zone. Because there is no refresh before
    // this method is called (at least during the very first cycle) - this can
    // including ghost zones can lead to timesteps of 0
    for (int iz = gz; iz < ngz - gz; iz++) {   // if rank < 3: gz = 0, ngz = 1
      for (int iy = gy; iy < ngy - gy; iy++) { // if rank < 2: gy = 0, ngy = 1
        for (int ix = gx; ix < ngx - gx; ix++) {
	  int i = INDEX(ix, iy, iz, ngx, ngy);
          dt = std::min(enzo_float(dt), std::abs(cooling_time[i]));
        }
      }
    }

    if (delete_cooling_time){
      delete [] cooling_time;
    }
  }
#endif

  return dt * courant_;
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE
void EnzoMethodGrackle::ResetEnergies ( EnzoBlock * enzo_block) throw()
{
   const EnzoConfig * enzo_config = enzo::config();
   EnzoUnits * enzo_units = enzo::units();

   chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

   /* Only need to do this if tracking chemistry */
   if (grackle_chemistry->primordial_chemistry < 1)
     return;

   Field field = enzo_block->data()->field();

   enzo_float * density     = (enzo_float*) field.values("density");
   enzo_float * internal_energy = (enzo_float*) field.values("internal_energy");
   enzo_float * total_energy    = (enzo_float*) field.values("total_energy");

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

   double temperature_slope = log10
     (enzo_config->initial_grackle_test_maximum_temperature/
      enzo_config->initial_grackle_test_minimum_temperature) / double(ny);

   for (int iz=gz; iz<nz+gz; iz++){ // Metallicity
     for (int iy=gy; iy<ny+gy; iy++) { // Temperature
       for (int ix=gx; ix<nx+gx; ix++) { // H Number Density
         int i = INDEX(ix,iy,iz,mx,my);

         enzo_float mu = e_density[i] + HI_density[i] + HII_density[i] +
            (HeI_density[i] + HeII_density[i] + HeIII_density[i])*0.25;

         if (grackle_chemistry->primordial_chemistry > 1){
           mu += HM_density[i] + 0.5 * (H2I_density[i] + H2II_density[i]);
         }

         if (grackle_chemistry->primordial_chemistry > 2){
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

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute_local_property_
(Block * block, enzo_float* values, code_units* grackle_units,
 grackle_field_data* grackle_fields, int i_hist,
 grackle_local_property_func func, std::string func_name) const throw()
{
  EnzoBlock * enzo_block = enzo::block(block);
  const EnzoConfig * enzo_config = enzo::config();

  code_units cur_grackle_units_;
  grackle_field_data cur_grackle_fields_;

  // setup grackle units if they are not already provided
  if (!grackle_units){
    grackle_units = &cur_grackle_units_;
    EnzoMethodGrackle::setup_grackle_units(enzo_block, grackle_units, i_hist);
  }

  // if grackle fields are not provided, define them
  bool delete_grackle_fields = false;
  if (!grackle_fields){
    grackle_fields  = &cur_grackle_fields_;
    EnzoMethodGrackle::setup_grackle_fields(enzo_block, grackle_fields,
					    i_hist);
    delete_grackle_fields = true;
  }

  // because this function is const-qualified, grackle_rates_ currently has
  // the type: const chemistry_data_storage *
  // we need to drop the `const` to be able to pass to to func (this is okay
  // because func will not actually modify the value).
  chemistry_data_storage * grackle_rates_ptr
    = const_cast<chemistry_data_storage *>(&grackle_rates_);

  if ((*func)(enzo_config->method_grackle_chemistry, grackle_rates_ptr,
	      grackle_units, grackle_fields, values) == ENZO_FAIL){
    ERROR1("EnzoMethodGrackle::compute_local_property_()",
	   "Error in call to Grackles's %s routine", func_name.c_str());
  }
  if (delete_grackle_fields){
    EnzoMethodGrackle::delete_grackle_fields(grackle_fields);
  }
  return;
}

#endif //CONFIG_USE_GRACKLE

//----------------------------------------------------------------------

void EnzoMethodGrackle::deallocate_grackle_rates_() throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  // sanity check:
  ASSERT("EnzoMethod::deallocate_grackle_rates_",
	 "enzo::config() must not return NULL",
	 enzo_config != NULL);
#ifdef CONFIG_USE_GRACKLE
  if (time_grackle_data_initialized_ == ENZO_FLOAT_UNDEFINED){
    ERROR("EnzoMethodGrackle::deallocate_grackle_rates_",
	  "grackle_rates_ data has not been allocated");
  }

  // deallocate previously the allocated allocated grackle_rates_ (doesn't
  // actually affect the chemistry_data pointer)
  _free_chemistry_data(enzo_config->method_grackle_chemistry, &grackle_rates_);
  // signal that grackle_data_ is not initialized
  time_grackle_data_initialized_ = ENZO_FLOAT_UNDEFINED;
#endif //CONFIG_USE_GRACKLE
}

//======================================================================
