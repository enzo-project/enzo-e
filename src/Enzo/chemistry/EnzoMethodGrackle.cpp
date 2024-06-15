// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
///           Andrew Emerick (aemerick11@gmail.com)
///           Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Tues May  7
/// @brief    Implements the EnzoMethodGrackle class

#include "Enzo/chemistry/chemistry.hpp"
#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"

//----------------------------------------------------------------------------

namespace { // anonymous namespace

// items within this namespace are just local to this file

GrackleChemistryData parse_chemistry(ParameterGroup p)
{
  // Initialize GrackleChemistryData
  // - we do this with a factory method that directly examines the parameter
  //   values within the "Method:grackle:*" group.
  // - Because Grackle has so many parameter values, it's very easy to make a
  //   small mistake when specifying the name of a parameter value and not
  //   notice until much later. For that reason, the factory method
  //   aggressively reports unexpected parameters as errors.
  // - to help this method, we provide 2 sets of parameter names

  //   1. specify all of the Grackle parameters that we will manually setup
  //      based on the values passed for other Enzo-E parameters. Errors will
  //      be reported if any of these are encountered
  const std::unordered_set<std::string> forbid_leaf_names = {"use_grackle",
                                                             "Gamma"};

  //   2. specify all parameters that MAY occur within the "Method:grackle:*"
  //      group that should be ignored by the factory method. (This needs to
  //      be updated if we introduce additional parameters for configuring
  //      EnzoMethodGrackle)
  const std::unordered_set<std::string> ignore_leaf_names =
    {"use_cooling_timestep", "radiation_redshift",
     // the next option is deprecated and is only listed in the short-term
     // for backwards compatability (it should now be replaced by
     // "Physics:fluid_props:floors:metallicity")
     "metallicity_floor",
     // for backwards compatability, we manually use "data_file" to
     // initialize "grackle_data_file" parameter (in the future, we may want
     // to change this)
     "data_file", "grackle_data_file",
     // the final two parameters auto-parsed by other Cello machinery
     "type", "courant"};

  GrackleChemistryData my_chemistry = GrackleChemistryData::from_parameters
    (p, forbid_leaf_names, ignore_leaf_names);

  // now let's manually initialize the handful of remaining runtime
  // parameters that are stored within method_grackle_chemistry

  // 1. use "Method:grackle:data_file" to initialize "grackle_data_file" 
  if (p.param("grackle_data_file") != nullptr){
    ERROR("EnzoMethodGrackle::from_parameters",
          "for backwards compatability, the user can't specify "
          "\"Method:grackle:grackle_data_file\". Instead, they must specify "
          "\"Method:grackle:data_file\".");
  } else if (p.param("data_file") != nullptr) {
    std::string fname = p.value_string("data_file", "");
    ASSERT("EnzoMethodGrackle::from_parameters",
           "\"Method:grackle:data_file\" can't be an empty string",
           fname.length() > 0); // sanity check!
    my_chemistry.set<std::string>("grackle_data_file", fname);
  } else {
    ERROR("EnzoMethodGrackle::from_parameters",
          "\"Method:grackle:data_file\" is required when using grackle");
  }

  // 2. update the value of use_grackle
  my_chemistry.set<int>("use_grackle", 1);

  // 3. Copy over parameters from Enzo-E to Grackle
  // (it's ok to check fluid-props since physics objects should be initialized
  //  before Method objects)
  if (enzo::fluid_props()->eos_variant().holds_alternative<EnzoEOSIdeal>()) {
    my_chemistry.set<double>
      ("Gamma", enzo::fluid_props()->eos_variant().get<EnzoEOSIdeal>().gamma);
  } else {
    ERROR("EnzoMethodGrackle::from_parameters",
          "Grackle currently can't be used when Enzo-E is configured to use "
          "an equation of state other than the ideal gas law");
  }

  // In the future, we may want to manually set use_radiative_transfer based
  // on an Enzo-E parameter for turning RT on / off:
  // my_chemistry.set<int>("use_radiative_transfer", ENZO_E_PARAMETER_NAME);

  return my_chemistry;
}

} // anonymous namespace

//----------------------------------------------------------------------------

EnzoMethodGrackle::EnzoMethodGrackle
(
 ParameterGroup p,
 const double physics_cosmology_initial_redshift,
 const double time
)
  : Method(),
    grackle_facade_(std::move(parse_chemistry(p)),
                    // for when not using cosmology - redshift of UVB
                    p.value_float("radiation_redshift", -1.0),
                    // the next parameter is relevant when using cosmology
                    physics_cosmology_initial_redshift,
                    time),
    use_cooling_timestep_(p.value_logical("use_cooling_timestep", false))
{
  // courant is only meaningful when use_cooling_timestep is true
  this->set_courant(p.value_float("courant", 1.0));

  // Gather list of fields that MUST be defined for this
  // method and check that they are permanent. If not,
  // define them.

  define_required_grackle_fields();

  /// Initialize default Refresh
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

}

//----------------------------------------------------------------------
 
void EnzoMethodGrackle::define_required_grackle_fields
()
{
  // Gather list of fields that MUST be defined for this method and
  // check that they are permanent. If not, define them.

  // This has been split off from the constructor so that other methods that
  // are initialized first and need knowledge of these fields at initialization
  // (e.g. to set up a refresh object), can ensure that the fields are defined

  const GrackleChemistryData* my_chemistry = this->try_get_chemistry();

  if ((my_chemistry == nullptr) ||
      (my_chemistry->get<int>("use_grackle") == 0)){
    return;
  }

  // special container for ensuring color fields are properly grouped

  std::vector<std::string> fields_to_define;
  std::vector<std::string> color_fields;

  const int rank = cello::rank();

  cello::define_field ("density");
  cello::define_field ("internal_energy");
  cello::define_field ("total_energy");

  if (rank>=1) cello::define_field ("velocity_x");
  if (rank>=2) cello::define_field ("velocity_y");
  if (rank>=3) cello::define_field ("velocity_z");

  // Get Grackle parameters defining fields to define

  const int metal_cooling = my_chemistry->get<int>("metal_cooling");
  const int chemistry_level = my_chemistry->get<int>("primordial_chemistry");
  const int specific_heating
    = my_chemistry->get<int>("use_specific_heating_rate");
  const int volumetric_heating
    = my_chemistry->get<int>("use_volumetric_heating_rate");
  const int radiative_transfer
      = my_chemistry->get<int>("use_radiative_transfer");

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
    cello::define_field_in_group ("DII_density", "color");
    cello::define_field_in_group ("HDI_density", "color");
  }

  if (radiative_transfer) {
    cello::define_field("RT_heating_rate");
    cello::define_field("RT_HI_ionization_rate");
    cello::define_field("RT_HeI_ionization_rate");
    cello::define_field("RT_HeII_ionization_rate");
    cello::define_field("RT_H2_dissociation_rate");
  }

  if (specific_heating) {
    cello::define_field("specific_heating_rate");
  }

  if (volumetric_heating) {
    cello::define_field("volumetric_heating_rate");
  }

}

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute ( Block * block) throw()
{

  if (block->is_leaf()){

#ifndef CONFIG_USE_GRACKLE

    ERROR("EnzoMethodGrackle::compute()",
          "Can't use method 'grackle' when Enzo-E isn't linked to Grackle");

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

//----------------------------------------------------------------------------

void EnzoMethodGrackle::update_grackle_density_fields(
                               Block * block,
                               grackle_field_data * grackle_fields
                               ) const throw() {
#ifndef CONFIG_USE_GRACKLE
  ERROR("EnzoMethodGrackle::update_grackle_density_fields",
        "Enzo-E isn't linked to grackle");
#else

  // Intended for use at problem initialization. Scale species
  // density fields to be sensible mass fractions of the initial
  // density field. Problem types that require finer-tuned control
  // over individual species fields should adapt this function
  // in their initialization routines.

  const GrackleChemistryData* my_chemistry = this->try_get_chemistry();
  ASSERT("update_grackle_density_fields", "not configured to use grackle",
         my_chemistry != nullptr);

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

  const int primordial_chemistry
    = my_chemistry->get<int>("primordial_chemistry");
  const int metal_cooling  = my_chemistry->get<int>("metal_cooling");
  const double HydrogenFractionByMass
    = my_chemistry->get<double>("HydrogenFractionByMass");
  const double DeuteriumToHydrogenRatio
    = my_chemistry->get<double>("DeuteriumToHydrogenRatio");

  for (int iz = 0; iz<ngz; iz++){
    for (int iy=0; iy<ngy; iy++){
      for (int ix=0; ix<ngx; ix++){
        int i = INDEX(ix,iy,iz,ngx,ngy);

        if(primordial_chemistry > 0){
          grackle_fields->HI_density[i]   = grackle_fields->density[i] * HydrogenFractionByMass;
          grackle_fields->HII_density[i]   = grackle_fields->density[i] * tiny_number;
          grackle_fields->HeI_density[i]   = grackle_fields->density[i] * (1.0 - HydrogenFractionByMass);
          grackle_fields->HeII_density[i]  = grackle_fields->density[i] * tiny_number;
          grackle_fields->HeIII_density[i] = grackle_fields->density[i] * tiny_number;
          grackle_fields->e_density[i]     = grackle_fields->density[i] * tiny_number;
        }

        if (primordial_chemistry > 1){
          grackle_fields->HM_density[i]    = grackle_fields->density[i] * tiny_number;
          grackle_fields->H2I_density[i]   = grackle_fields->density[i] * tiny_number;
          grackle_fields->H2II_density[i]  = grackle_fields->density[i] * tiny_number;
        }

        if (primordial_chemistry > 2){
          grackle_fields->DI_density[i]    = grackle_fields->density[i] * DeuteriumToHydrogenRatio;
          grackle_fields->DII_density[i]   = grackle_fields->density[i] * tiny_number;
          grackle_fields->HDI_density[i]   = grackle_fields->density[i] * tiny_number;
        }

       if (metal_cooling == 1){
          grackle_fields->metal_density[i] = grackle_fields->density[i] * metal_factor;
       }

      }
    }
  }

  if (cleanup_grackle_fields){
    EnzoMethodGrackle::delete_grackle_fields(grackle_fields);
  }

  return;

#endif // CONFIG_USE_GRACKLE
}

//----------------------------------------------------------------------

void EnzoMethodGrackle::compute_ ( Block * block) throw()
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("EnzoMethodGrackle::compute_", "Enzo-E isn't linked to grackle");
#else
  const EnzoConfig * enzo_config = enzo::config();
  if (cello::is_initial_cycle(InitCycleKind::fresh_or_noncharm_restart)) {
    bool nohydro = ( (enzo::problem()->method("ppm") == nullptr) |
                     (enzo::problem()->method("mhd_vlct") == nullptr) |
                     (enzo::problem()->method("ppml") == nullptr) );

    ASSERT("EnzoMethodGrackle::compute_",
           "The current implementation requires the dual-energy formalism to "
           "be in use, when EnzoMethodGrackle is used with a (M)HD-solver",
           nohydro | !enzo::fluid_props()->dual_energy_config().is_disabled());
  }

  // Solve chemistry
  // NOTE: should we set compute_time to `time + 0.5*dt`?
  //       I think that's what enzo-classic does...
  double compute_time = block->state()->time(); // only matters in cosmological sims
  grackle_facade_.solve_chemistry(block, compute_time, block->state()->dt());

  // now we have to do some extra-work after the fact (such as adjusting total
  // energy density and applying floors...)

  // todo: avoid constructing this instance of grackle_fields
  grackle_field_data grackle_fields;
  setup_grackle_fields(EnzoFieldAdaptor(block,0), &grackle_fields);

  Field field = block->data()->field();

  int gx,gy,gz;
  field.ghost_depth (0,&gx,&gy,&gz);

  int nx,ny,nz;
  field.size (&nx,&ny,&nz);

  int ngx = nx + 2*gx;
  int ngy = ny + 2*gy;
  int ngz = nz + 2*gz;

  const int rank = cello::rank();

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

  delete_grackle_fields(&grackle_fields);

  return;
#endif // CONFIG_USE_GRACKLE
}

//----------------------------------------------------------------------

double EnzoMethodGrackle::timestep ( Block * block ) throw()
{
  double dt = std::numeric_limits<double>::max();;

  if (use_cooling_timestep_){
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

    this->calculate_cooling_time(EnzoFieldAdaptor(block,0), cooling_time, 0);

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

  return dt * courant_;
}

//----------------------------------------------------------------------

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
