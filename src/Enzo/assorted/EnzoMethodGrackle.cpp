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
    ,
    grackle_units_(),
    grackle_rates_(),
    time_grackle_data_initialized_(ENZO_FLOAT_UNDEFINED)
#endif
{
#ifdef CONFIG_USE_GRACKLE

  // Gather list of fields that MUST be defined for this
  // method and check that they are permanent. If not,
  // define them.

  define_required_grackle_fields();

  /// Initialize default Refresh
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  /// Define Grackle's internal data structures
  time_grackle_data_initialized_ = ENZO_FLOAT_UNDEFINED;
  initialize_grackle_chemistry_data(time);
#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE
 
void EnzoMethodGrackle::define_required_grackle_fields()
{
  // Gather list of fields that MUST be defined for this method and
  // check that they are permanent. If not, define them.

  // This has been split off from the constructor so that other methods that
  // are initialized first and need knowledge of these fields at initialization
  // (e.g. to set up a refresh object), can ensure that the fields are defined

  if (!enzo::config()->method_grackle_use_grackle) {return;}

  // special container for ensuring color fields are properly grouped

  std::vector<std::string> fields_to_define;
  std::vector<std::string> color_fields;

  chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

  const int rank = cello::rank();

  cello::define_field ("density");
  cello::define_field ("internal_energy");
  cello::define_field ("total_energy");

  if (rank>=1) cello::define_field ("velocity_x");
  if (rank>=2) cello::define_field ("velocity_y");
  if (rank>=3) cello::define_field ("velocity_z");

  // Get Grackle parameters defining fields to define

  const int metal_cooling      = grackle_chemistry->metal_cooling;
  const int chemistry_level    = grackle_chemistry->primordial_chemistry;
  const int specific_heating   = grackle_chemistry->use_specific_heating_rate;
  const int volumetric_heating = grackle_chemistry->use_volumetric_heating_rate;

  // Metal cooling fields

  if (metal_cooling > 0) {
    cello::define_field_in_group ("metal_density", "color");
  }

  // Primordial chemistry fields

  if (chemistry_level >= 1) {
    cello::define_field_in_group ("HI_density",    "color");
    cello::define_field_in_group ("HII_density",   "color");
    cello::define_field_in_group ("HeI_density",   "color");
    cello::define_field_in_group ("HeII_density",  "color");
    cello::define_field_in_group ("HeIII_density", "color");
    cello::define_field_in_group ("e_density",     "color");
  }

  if (chemistry_level >= 2) {
    cello::define_field_in_group ("HM_density",   "color");
    cello::define_field_in_group ("H2I_density",  "color");
    cello::define_field_in_group ("H2II_density", "color");
  }

  if (chemistry_level >= 3) {
    cello::define_field_in_group ("DI_density",  "color");
    cello::define_field_in_group ("DII_density", "color" );
    cello::define_field_in_group ("HDI_density", "color");
  }

  if (specific_heating) {
    cello::define_field("specific_heating_rate");
  }

  if (volumetric_heating) {
    cello::define_field("volumetric_heating_rate");
  }

}
#endif /* CONFIG_USE_GRACKLE */

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute ( Block * block) throw()
{

  if (block->is_leaf()){

#ifndef CONFIG_USE_GRACKLE

    ERROR("EnzoMethodGrackle::compute()",
          "Trying to use method 'grackle' with "
          "Grackle configuration turned off!");

#else /* CONFIG_USE_GRACKLE */

    // Start timer
    Simulation * simulation = cello::simulation();
    if (simulation)
      simulation->performance()->start_region(perf_grackle,__FILE__,__LINE__);

    this->compute_(block);

    if (simulation)
      simulation->performance()->stop_region(perf_grackle,__FILE__,__LINE__);
#endif
  }

  block->compute_done();

  return;

}

//----------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE

void EnzoMethodGrackle::initialize_grackle_chemistry_data
(double current_time, bool preinitialized_units /* default false */)
{

  /* Define Grackle's chemistry data if not yet defined */

  if (this->time_grackle_data_initialized_ == current_time) return;

  if (this->time_grackle_data_initialized_ != ENZO_FLOAT_UNDEFINED){
    // deallocate previously the allocated allocated grackle_rates_ (doesn't
    // actually affect the chemistry_data pointer)
    deallocate_grackle_rates_();
  }

  // Initialize grackle units and data
  TRACE("Calling initialize_chemistry_data from EnzoMethodGrackle::EnzoMethodGrackle() or EnzoMethodGrackle::pup(PUP::er &p)");

  const EnzoConfig * enzo_config = enzo::config();

  if (!preinitialized_units){
    // initialize grackle_units
    setup_grackle_units (current_time, &grackle_units_);
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

#endif

//----------------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE

void EnzoMethodGrackle::setup_grackle_units (double current_time,
                                             code_units * grackle_units
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
    if (current_time < 0){
      ERROR("EnzoMethodGrackle::setup_grackle_units",
            "A valid current_time value is required");
    }

    enzo_float cosmo_a  = 1.0;
    enzo_float cosmo_dt = 0.0;

    EnzoPhysicsCosmology * cosmology = enzo::cosmology();

    grackle_units->density_units  = cosmology->density_units();
    grackle_units->length_units = cosmology->length_units();
    grackle_units->time_units  = cosmology->time_units();
    grackle_units->velocity_units = cosmology->velocity_units();

    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dt, current_time);
    grackle_units->a_units
         = 1.0 / (1.0 + enzo_config->physics_cosmology_initial_redshift);
    grackle_units->a_value = cosmo_a;

  } else if (enzo_config->method_grackle_radiation_redshift > -1){
    grackle_units->a_value = 1.0 /
                         (1.0 + enzo_config->method_grackle_radiation_redshift);
  }

  return;
}

#endif

//----------------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE

void EnzoMethodGrackle::setup_grackle_units (const EnzoFieldAdaptor& fadaptor,
                                             code_units * grackle_units
                                             ) throw()
{
  const EnzoConfig * config = enzo::config();
    double current_time =
      (config->physics_cosmology) ? fadaptor.compute_time() : -1.0;
    setup_grackle_units(current_time, grackle_units);
}

#endif

//--------------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE

void EnzoMethodGrackle::setup_grackle_fields
(const EnzoFieldAdaptor& fadaptor,
 grackle_field_data * grackle_fields,
 int stale_depth, /* default: 0 */
 bool omit_cell_width /* default false */
 ) const throw()
{

  // Grackle grid dimenstion and grid size
  grackle_fields->grid_rank      = cello::rank();
  grackle_fields->grid_dimension = new int[3];
  grackle_fields->grid_start     = new int[3];
  grackle_fields->grid_end       = new int[3];

  fadaptor.get_grackle_field_grid_props(grackle_fields->grid_dimension,
                                        grackle_fields->grid_start,
                                        grackle_fields->grid_end);

  if (stale_depth > 0){
    ERROR("EnzoMethodGrackle::setup_grackle_fields", "untested");
    for (int i = 1; i <= grackle_fields->grid_rank; i++){
      grackle_fields->grid_start[i-1] += stale_depth;
      grackle_fields->grid_end[i-1] -= stale_depth;

      // reminder for following check: grackle_fields->grid_end is inclusive
      if (grackle_fields->grid_end[i-1] < grackle_fields->grid_start[i-1]){
        ERROR("EnzoMethodGrackle::setup_grackle_fields",
              "stale_depth is too large");
      }
    }
  } else if (stale_depth < 0){
    ERROR("EnzoMethodGrackle::setup_grackle_fields",
          "can't handle negative stale_depth");
  }

  if (omit_cell_width){
    grackle_fields->grid_dx = 0.0;
  } else {
    double hx, hy, hz;
    fadaptor.cell_width(&hx,&hy,&hz);
    grackle_fields->grid_dx = hx;
  }

  // Setup all fields to be passed into grackle
  grackle_fields->density         = fadaptor.ptr_for_grackle("density", true);
  grackle_fields->internal_energy = fadaptor.ptr_for_grackle("internal_energy",
                                                             true);
  grackle_fields->x_velocity      = fadaptor.ptr_for_grackle("velocity_x");
  grackle_fields->y_velocity      = fadaptor.ptr_for_grackle("velocity_y");
  grackle_fields->z_velocity      = fadaptor.ptr_for_grackle("velocity_z");

  // Get chemical species fields if they exist

  // primordial_chemistry > 0 fields
  grackle_fields->HI_density      = fadaptor.ptr_for_grackle("HI_density");
  grackle_fields->HII_density     = fadaptor.ptr_for_grackle("HII_density");
  grackle_fields->HeI_density     = fadaptor.ptr_for_grackle("HeI_density");
  grackle_fields->HeII_density    = fadaptor.ptr_for_grackle("HeII_density");
  grackle_fields->HeIII_density   = fadaptor.ptr_for_grackle("HeIII_density");
  grackle_fields->e_density       = fadaptor.ptr_for_grackle("e_density");

  // primordial_chemistry > 1 fields
  grackle_fields->HM_density      = fadaptor.ptr_for_grackle("HM_density");
  grackle_fields->H2I_density     = fadaptor.ptr_for_grackle("H2I_density");
  grackle_fields->H2II_density    = fadaptor.ptr_for_grackle("H2II_density");

  // primordial_chemistry > 2 fields
  grackle_fields->DI_density      = fadaptor.ptr_for_grackle("DI_density");
  grackle_fields->DII_density     = fadaptor.ptr_for_grackle("DII_density");
  grackle_fields->HDI_density     = fadaptor.ptr_for_grackle("HDI_density");

  grackle_fields->metal_density   = fadaptor.ptr_for_grackle("metal_density");

  /* Leave these as NULL for now and save for future development */
  gr_float * volumetric_heating_rate = NULL;
  gr_float * specific_heating_rate = NULL;

  grackle_fields->volumetric_heating_rate = volumetric_heating_rate;
  grackle_fields->specific_heating_rate   = specific_heating_rate;

  return;
}

#endif

//----------------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE

void EnzoMethodGrackle::update_grackle_density_fields(
                               Block * block,
                               grackle_field_data * grackle_fields
                               ) const throw() {

  // Intended for use at problem initialization. Scale species
  // density fields to be sensible mass fractions of the initial
  // density field. Problem types that require finer-tuned control
  // over individual species fields should adapt this function
  // in their initialization routines.

  grackle_field_data tmp_grackle_fields;
  bool cleanup_grackle_fields = false;
  if (grackle_fields == nullptr){
    grackle_fields = &tmp_grackle_fields;
    setup_grackle_fields(block, grackle_fields);
    cleanup_grackle_fields = true;
  }

  Field field = block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  const double tiny_number = 1.0E-10;

  const EnzoFluidFloorConfig& fluid_floors
    = enzo::fluid_props()->fluid_floor_config();
  const enzo_float metal_factor = fluid_floors.has_metal_mass_frac_floor()
    ? fluid_floors.metal_mass_frac() : (enzo_float)tiny_number;

  const EnzoConfig * enzo_config = enzo::config();
  chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;
  double metallicity_floor_ =
    enzo::fluid_props()->fluid_floor_config().metal_mass_frac();

  for (int iz = 0; iz<ngz; iz++){
    for (int iy=0; iy<ngy; iy++){
      for (int ix=0; ix<ngx; ix++){
        int i = INDEX(ix,iy,iz,ngx,ngy);

        if(grackle_chemistry->primordial_chemistry > 0){
          grackle_fields->HI_density[i]   = grackle_fields->density[i] * grackle_chemistry->HydrogenFractionByMass;
          grackle_fields->HII_density[i]   = grackle_fields->density[i] * tiny_number;
          grackle_fields->HeI_density[i]   = grackle_fields->density[i] * (1.0 - grackle_chemistry->HydrogenFractionByMass);
          grackle_fields->HeII_density[i]  = grackle_fields->density[i] * tiny_number;
          grackle_fields->HeIII_density[i] = grackle_fields->density[i] * tiny_number;
          grackle_fields->e_density[i]     = grackle_fields->density[i] * tiny_number;
        }

        if (grackle_chemistry->primordial_chemistry > 1){
          grackle_fields->HM_density[i]    = grackle_fields->density[i] * tiny_number;
          grackle_fields->H2I_density[i]   = grackle_fields->density[i] * tiny_number;
          grackle_fields->H2II_density[i]  = grackle_fields->density[i] * tiny_number;
        }

        if (grackle_chemistry->primordial_chemistry > 2){
          grackle_fields->DI_density[i]    = grackle_fields->density[i] * grackle_chemistry->DeuteriumToHydrogenRatio;
          grackle_fields->DII_density[i]   = grackle_fields->density[i] * tiny_number;
          grackle_fields->HDI_density[i]   = grackle_fields->density[i] * tiny_number;
        }

       if (grackle_chemistry->metal_cooling == 1){
          grackle_fields->metal_density[i] = grackle_fields->density[i] * metal_factor;
       }

      }
    }
  }

  if (cleanup_grackle_fields){
    EnzoMethodGrackle::delete_grackle_fields(grackle_fields);
  }

  return;
}

#endif

//----------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE

void EnzoMethodGrackle::compute_ ( Block * block) throw()
{
  const EnzoConfig * enzo_config = enzo::config();
  if (block->cycle() == enzo_config->initial_cycle) {
    bool nohydro = ( (enzo::problem()->method("ppm") == nullptr) |
                     (enzo::problem()->method("mhd_vlct") == nullptr) |
                     (enzo::problem()->method("ppml") == nullptr) );

    ASSERT("EnzoMethodGrackle::compute_",
           "The current implementation requires the dual-energy formalism to "
           "be in use, when EnzoMethodGrackle is used with a (M)HD-solver",
           nohydro | !enzo::fluid_props()->dual_energy_config().is_disabled());
  }

  Field field = block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  const int rank = cello::rank();

  /* Set code units for use in grackle */
  grackle_field_data grackle_fields;

  EnzoFieldAdaptor fadaptor(block, 0);

  setup_grackle_units(fadaptor, &this->grackle_units_);
  setup_grackle_fields(fadaptor, &grackle_fields);

  chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

  // Solve chemistry
  double dt = block->dt();
  if (local_solve_chemistry(grackle_chemistry, &grackle_rates_,
			    &grackle_units_, &grackle_fields, dt)
      == ENZO_FAIL) {
    ERROR("EnzoMethodGrackle::compute()",
    "Error in local_solve_chemistry.\n");
  }

  // enforce metallicity floor (if one was provided)
  enforce_metallicity_floor(block);

  /* Correct total energy for changes in internal energy */
  gr_float * v3[3];
  v3[0] = grackle_fields.x_velocity;
  v3[1] = grackle_fields.y_velocity;
  v3[2] = grackle_fields.z_velocity;

  const bool mhd = field.is_field("bfield_x");
  enzo_float * b3[3] = {NULL, NULL, NULL};
  if (mhd) {
    b3[0]                = (enzo_float*) field.values("bfield_x");
    if (rank >= 2) b3[1] = (enzo_float*) field.values("bfield_y");
    if (rank >= 3) b3[2] = (enzo_float*) field.values("bfield_z");
  }

  enzo_float * total_energy    = (enzo_float *) field.values("total_energy");
  for (int i = 0; i < ngx*ngy*ngz; i++){
    total_energy[i] = grackle_fields.internal_energy[i];

    enzo_float inv_density;
    if (mhd) inv_density = 1.0 / grackle_fields.density[i];
    for (int dim = 0; dim < rank; dim++){
      total_energy[i] += 0.5 * v3[dim][i] * v3[dim][i];
      if (mhd) total_energy[i] += 0.5 * b3[dim][i] * b3[dim][i] * inv_density;
    }
  }

  // For testing purposes - reset internal energies with changes in mu
  if (enzo_config->initial_grackle_test_reset_energies){
    this->ResetEnergies(block);
  }

  delete_grackle_fields(&grackle_fields);

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
    Field field = block->data()->field();

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

    calculate_cooling_time(EnzoFieldAdaptor(block,0), cooling_time, 0,
                           nullptr, nullptr);

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

void EnzoMethodGrackle::enforce_metallicity_floor(Block * block) throw()
{
  const EnzoFluidFloorConfig& fluid_floors
    = enzo::fluid_props()->fluid_floor_config();

  if (!fluid_floors.has_metal_mass_frac_floor()){
    return; // return early if the floor has not been defined
  }

  const enzo_float metal_mass_frac_floor = fluid_floors.metal_mass_frac();

  // MUST have metal_density field tracked
  Field field = block->data()->field();
  enzo_float * density = (enzo_float*) field.values("density");
  enzo_float * metal_density  = (enzo_float*) field.values("metal_density");
  ASSERT("EnzoMethodGrackle::enforce_metallicity_floor",
         ("Can't enforce metallicity floor when the \"metal_density\" field "
          "doesn't exist"),
         metal_density != nullptr);

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  for (int iz=0; iz<ngz; iz++){
    for (int iy=0; iy<ngy; iy++){
      for (int ix=0; ix<ngx; ix++){
        int i = INDEX(ix,iy,iz,ngx,ngy);
        metal_density[i] = std::max(metal_density[i],
                                    metal_mass_frac_floor * density[i]);
      }
    }
  }
  return;
}
#endif
//----------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE
void EnzoMethodGrackle::ResetEnergies ( Block * block) throw()
{
   const EnzoConfig * enzo_config = enzo::config();
   EnzoUnits * enzo_units = enzo::units();

   chemistry_data * grackle_chemistry =
    enzo::config()->method_grackle_chemistry;

   /* Only need to do this if tracking chemistry */
   if (grackle_chemistry->primordial_chemistry < 1)
     return;

   Field field = block->data()->field();

   enzo_float * density     = (enzo_float*) field.values("density");
   enzo_float * internal_energy = (enzo_float*) field.values("internal_energy");
   enzo_float * total_energy    = (enzo_float*) field.values("total_energy");

   enzo_float * HI_density    = (enzo_float*) field.values("HI_density");
   enzo_float * HII_density   = (enzo_float*) field.values("HII_density");
   enzo_float * HeI_density   = (enzo_float*) field.values("HeI_density");
   enzo_float * HeII_density  = (enzo_float*) field.values("HeII_density");
   enzo_float * HeIII_density = (enzo_float*) field.values("HeIII_density");
   enzo_float * e_density     = (enzo_float*) field.values("e_density");

   enzo_float * H2I_density   = (enzo_float*) field.values("H2I_density");
   enzo_float * H2II_density  = (enzo_float*) field.values("H2II_density");
   enzo_float * HM_density    = (enzo_float*) field.values("HM_density");

   enzo_float * DI_density    = (enzo_float *) field.values("DI_density");
   enzo_float * DII_density   = (enzo_float *) field.values("DII_density");
   enzo_float * HDI_density   = (enzo_float *) field.values("HDI_density");

   // JB: NOTE note metal_density field accessed
   enzo_float * metal_density = (enzo_float*) field.values("metal_density");

   // Field size
   int nx,ny,nz;
   field.size(&nx,&ny,&nz);

   // Cell widths
   double xm,ym,zm;
   block->data()->lower(&xm,&ym,&zm);
   double xp,yp,zp;
   block->data()->upper(&xp,&yp,&zp);

   // Ghost depths
   int gx,gy,gz;
   field.ghost_depth(0,&gx,&gy,&gz);

   int mx,my,mz;
   field.dimensions(0,&mx,&my,&mz);

   double temperature_slope = log10
     (enzo_config->initial_grackle_test_maximum_temperature/
      enzo_config->initial_grackle_test_minimum_temperature) / double(ny);

   const enzo_float nominal_gamma = enzo::fluid_props()->gamma();

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

         internal_energy[i] =
           (pow(10.0, ((temperature_slope * (iy-gy)) +
                       log10(enzo_config->initial_grackle_test_minimum_temperature)))/
            mu / enzo_units->kelvin_per_energy_units() /
            (nominal_gamma - 1.0));
         total_energy[i] = internal_energy[i];

       }
     }
   }

  return;
}

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute_local_property_
(const EnzoFieldAdaptor& fadaptor, enzo_float* values, int stale_depth,
 code_units* grackle_units, grackle_field_data* grackle_fields,
 grackle_local_property_func func, std::string func_name) const throw()
{
  const EnzoConfig * enzo_config = enzo::config();

  code_units cur_grackle_units_;
  grackle_field_data cur_grackle_fields;

  // setup grackle units if they are not already provided
  if (!grackle_units){
    grackle_units = &cur_grackle_units_;
    EnzoMethodGrackle::setup_grackle_units(fadaptor, grackle_units);
  }

  // if grackle fields are not provided, define them
  bool delete_grackle_fields = false;
  if (!grackle_fields){
    // the cell width is not used for computing for local properties. Thus, we
    // don't require it (for cases where fadaptor wraps EnzoEFltArrayMap)
    bool omit_cell_width = true;

    grackle_fields  = &cur_grackle_fields;
    EnzoMethodGrackle::setup_grackle_fields(fadaptor, grackle_fields,
                                            stale_depth, omit_cell_width);
    delete_grackle_fields = true;
  }

  for (int i = 1; i <= grackle_fields->grid_rank; i++){
    int ax_start = grackle_fields->grid_start[i-1];
    int ax_end = grackle_fields->grid_end[i-1];
    int ax_dim = grackle_fields->grid_dimension[i-1];

    // currently, Grackle's local_calculate_pressure, local_calculate_gamma,
    // & local_calculate_temperature functions ignore grid_start & grid_end
    // This assertion makes sure users won't get unexpected results...
    ASSERT("EnzoMethodGrackle::compute_local_property_",
           ("until PR #106 is merged into Grackle, we require "
            "grackle_fields->grid_start & grackle_fields->grid_end to include "
            "all data"),
           (ax_start == 0) & ((ax_end+1) == ax_dim));
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
