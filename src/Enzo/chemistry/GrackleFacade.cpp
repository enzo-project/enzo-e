// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_GrackleFacade.cpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Feb 27 2023
/// @brief    [\ref Enzo] Implementation of GrackleFacade class
///
/// This class provides a facade to the Grackle library. It is responsible for
/// managing the lifetime of all associated objects

#include "cello.hpp"
#include "enzo.hpp"

//----------------------------------------------------------------------------

#ifndef CONFIG_USE_GRACKLE
// provide dummy definitions of code_units and chemistry_data_storage
// - these are both required by std::unique_ptr
// - the dummy definition of code_units also comes in handy for defining some
//   helper functions that are local to this file
extern "C" {
  struct code_units { int dummy; };
  struct chemistry_data_storage { int dummy; };
}
#endif

//----------------------------------------------------------------------------

bool GrackleFacade::linked_against_grackle() noexcept
{
#ifdef CONFIG_USE_GRACKLE
  return true;
#else
  return false;
#endif
}

//----------------------------------------------------------------------------

namespace { // things within anonymous namespace are local to this file

// PUP operator for Grackle's code_units
void operator|(PUP::er &p, code_units *c){
#ifndef CONFIG_USE_GRACKLE
  ERROR("operator|",
        "can't define PUP operation for code_units when Grackle isn't used");
#else

  ASSERT("operator|", "ptr to code_units struct can't be a nullptr",
         c != nullptr);

  // Make sure to change this if code_units ever changes
  // all are just single values (int or float)

  p | c->comoving_coordinates;
  p | c->density_units;
  p | c->length_units;
  p | c->time_units;
  p | c->velocity_units;
  p | c->a_units;
  p | c->a_value;

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------------

/// returns a pointer to a fully setup a code_units struct
///
/// When grackle_units is a nullptr, this allocates the returned code_units
/// object on the heap. Otherwise, the configured/returned object is the one
/// that was passed to the grackle_units argument
///
/// It's worth discussing a weird quirk in Grackle's API when it comes to the
/// ``code_units`` struct:
///  - as a general rule, after initializing Grackle none of ``code_units``'s
///    fields should EVER be mutated (if they change, there's a strong chance
///    that the units data probably becomes invalidated).
///  - The only exception to this rule is the ``a_value`` field when using
///    comoving coordinates (this affects the quantities like the UV
///    background)
/// In principle, one could mutate the initial ``code_units`` struct or mutate
/// a copy of the original. However, we currently just rebuild it from scratch
/// every time we need it.
code_units* setup_grackle_u_(double current_time, double radiation_redshift,
                             code_units* grackle_units = nullptr) noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("setup_grackle_units", "grackle isn't being used");
#else

  if (grackle_units == nullptr) {
    grackle_units = new code_units;
  }

  EnzoUnits * enzo_units = enzo::units();
  EnzoPhysicsCosmology * cosmology = enzo::cosmology();

  /* Set code units for use in grackle */

  grackle_units->comoving_coordinates = (cosmology != nullptr);
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

    // Historically, there have been some issues with enzo_units being properly
    // configured to know that it is a cosmological simulation by time this
    // function is called. For that reason, we take the units directly from the
    // cosmology object.
    grackle_units->density_units  = cosmology->density_units();
    grackle_units->length_units = cosmology->length_units();
    grackle_units->time_units  = cosmology->time_units();
    grackle_units->velocity_units = cosmology->velocity_units();

    cosmology->compute_expansion_factor(&cosmo_a, &cosmo_dt, current_time);
    grackle_units->a_units = 1.0 / (1.0 + cosmology->initial_redshift());
    grackle_units->a_value = cosmo_a;

  } else if (radiation_redshift >= 0){
    grackle_units->a_value = 1.0 / (1.0 + radiation_redshift);
  }

  return grackle_units;
#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------------

std::unique_ptr<chemistry_data_storage> init_rates_
(GrackleChemistryData& my_chemistry, code_units* grackle_units)
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("init_rates_", "Grackle isn't linked");
#else

  std::unique_ptr<chemistry_data_storage> out(new chemistry_data_storage);

  if (local_initialize_chemistry_data(my_chemistry.get_ptr(),
                                      out.get(), grackle_units) == ENZO_FAIL) {
    ERROR("init_rates_", "Error in _initialize_chemistry_data");
  }

  return out;
#endif /* CONFIG_USE_GRACKLE */
}

} // anonymous namespace

//----------------------------------------------------------------------------

GrackleFacade::GrackleFacade(GrackleChemistryData&& my_chemistry,
                             const double radiation_redshift,
                             const double physics_cosmology_initial_redshift,
                             const double time) throw()
  : my_chemistry_(std::move(my_chemistry)),
    grackle_units_(nullptr),
    grackle_rates_(nullptr),
    radiation_redshift_(radiation_redshift)
{
  if ((radiation_redshift >= 0) && (enzo::cosmology() != nullptr)){
    ERROR("GrackleFacade::GrackleFacade",
          "radiation redshift should not be specified in cosmological sims");
  }
  grackle_units_ = std::unique_ptr<code_units>
    (setup_grackle_u_(time,radiation_redshift, nullptr));
  grackle_rates_ = init_rates_(my_chemistry_, grackle_units_.get());
}

//----------------------------------------------------------------------------

GrackleFacade::GrackleFacade(CkMigrateMessage *m)
  : PUP::able(m),
    my_chemistry_(),
    grackle_units_(nullptr),
    grackle_rates_(nullptr),
    radiation_redshift_(-1)
{ }

//----------------------------------------------------------------------------

void GrackleFacade::pup(PUP::er &p) {
  const bool up = p.isUnpacking();
  if (up && ((grackle_units_.get() != nullptr) ||
             (grackle_rates_.get() != nullptr)) ){
    ERROR("GrackleFacade::pup",
          "can't unpack data into an already initialized instance");
  }

  // NOTE: initializing grackle_units_ from scratch requires that
  //       enzo::physics() returns a properly configured global physics object.
  //       We choose to pup grackle_units_ directly so that we don't need to
  //       worry about whether that global object has already been unpacked.

  PUP::able::pup(p);
  p | my_chemistry_;

  if (up) {
    this->grackle_units_ = std::unique_ptr<code_units>(new code_units);
  }
  p | grackle_units_.get();
  p | radiation_redshift_;

  if (up) {
    grackle_rates_ = init_rates_(my_chemistry_, grackle_units_.get());
  }
}

//----------------------------------------------------------------------------

GrackleFacade::~GrackleFacade() noexcept
{
  // deallocate the previously allocated fields of grackle_rates_. This doesn't
  // do either of the following:
  // - affect the contents of my_chemistry_ (that's handled later in its
  //   destructor)
  // - deallocate the grackle_rates_ pointer (that's handled later by the
  //   destructor of std::unique_ptr)
#ifdef CONFIG_USE_GRACKLE

  // the ONLY way grackle_rates could be a nullptr is if something went wrong
  // and the pup method never got called after the migration constructor
  // (in that case the program would already have aborted)
  local_free_chemistry_data(my_chemistry_.get_ptr(), grackle_rates_.get());

#endif
}

//----------------------------------------------------------------------------

void GrackleFacade::setup_grackle_fields
(const EnzoFieldAdaptor& fadaptor,
 grackle_field_data * grackle_fields,
 int stale_depth, /* default: 0 */
 bool omit_cell_width /* default false */
 ) const noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleFacade::setup_grackle_fields", "grackle isn't being used");
#else

  // Grackle grid dimenstion and grid size
  grackle_fields->grid_rank      = cello::rank();
  grackle_fields->grid_dimension = new int[3];
  grackle_fields->grid_start     = new int[3];
  grackle_fields->grid_end       = new int[3];

  fadaptor.get_grackle_field_grid_props(grackle_fields->grid_dimension,
                                        grackle_fields->grid_start,
                                        grackle_fields->grid_end);

  if (stale_depth > 0){
    ERROR("GrackleFacade::setup_grackle_fields", "untested");
    for (int i = 1; i <= grackle_fields->grid_rank; i++){
      grackle_fields->grid_start[i-1] += stale_depth;
      grackle_fields->grid_end[i-1] -= stale_depth;

      // reminder for following check: grackle_fields->grid_end is inclusive
      if (grackle_fields->grid_end[i-1] < grackle_fields->grid_start[i-1]){
        ERROR("GrackleFacade::setup_grackle_fields",
              "stale_depth is too large");
      }
    }
  } else if (stale_depth < 0){
    ERROR("GrackleFacade::setup_grackle_fields",
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

  //  RADIATIVE TRANSFER HEATING AND IONIZATION RATES
  grackle_fields->RT_heating_rate         = fadaptor.ptr_for_grackle("RT_heating_rate");
  grackle_fields->RT_HI_ionization_rate   = fadaptor.ptr_for_grackle("RT_HI_ionization_rate");
  grackle_fields->RT_HeI_ionization_rate  = fadaptor.ptr_for_grackle("RT_HeI_ionization_rate");
  grackle_fields->RT_HeII_ionization_rate = fadaptor.ptr_for_grackle("RT_HeII_ionization_rate");
  grackle_fields->RT_H2_dissociation_rate = fadaptor.ptr_for_grackle("RT_H2_dissociation_rate");

  /* Leave these as NULL for now and save for future development */
  gr_float * volumetric_heating_rate = NULL;
  gr_float * specific_heating_rate = NULL;

  grackle_fields->volumetric_heating_rate = volumetric_heating_rate;
  grackle_fields->specific_heating_rate   = specific_heating_rate;

#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------------

void GrackleFacade::delete_grackle_fields(grackle_field_data * grackle_fields)
  const noexcept
{
#ifdef CONFIG_USE_GRACKLE
  // previously versions set members to nullptr, but that's not strictly
  // necessary AND it was non-exhaustive
  delete [] grackle_fields->grid_dimension;
  delete [] grackle_fields->grid_start;
  delete [] grackle_fields->grid_end;
#endif /* CONFIG_USE_GRACKLE */
}

//----------------------------------------------------------------------------

void GrackleFacade::solve_chemistry(Block* block, double compute_time,
                                    double dt) const noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleFacade::solve_chemistry", "grackle isn't being used");
#else

  code_units grackle_units;
  setup_grackle_u_(compute_time, radiation_redshift_, &grackle_units);

  EnzoFieldAdaptor fadaptor(block, 0);
  grackle_field_data grackle_fields;
  setup_grackle_fields(fadaptor, &grackle_fields);

  // because this function is const-qualified, my_chemistry_.get_ptr()
  // currently returns a pointer to a `const`. we need to drop the `const` to
  // be able to pass it to local_solve_chemistry (this is ok because the
  // function won't actually modify the value).
  chemistry_data * chemistry_data_ptr
    = const_cast<chemistry_data *>(my_chemistry_.get_ptr());

  if (local_solve_chemistry(chemistry_data_ptr, grackle_rates_.get(),
			    &grackle_units, &grackle_fields, dt)
      == ENZO_FAIL) {
    ERROR("GrackleFacade::solve_chemistry", "Error in local_solve_chemistry.");
  }

  delete_grackle_fields(&grackle_fields);
#endif
}

//----------------------------------------------------------------------------

typedef int (*grackle_local_property_func)(chemistry_data*,
                                           chemistry_data_storage *,
                                           code_units*, grackle_field_data*,
                                           enzo_float*);

//----------------------------------------------------------------------------

namespace { // things within anonymous namespace are local to this file

struct GracklePropFnInfo{
  grackle_local_property_func ptr;
  const char* name_suffix;
  bool has_z_dependence;
};

//----------------------------------------------------------------------------

#ifdef CONFIG_USE_GRACKLE
GracklePropFnInfo get_fn_info_(GracklePropertyEnum func_choice)
{
  // the value of the has_z_dependence field was determined by inspecting
  // Grackle's source code.
  switch (func_choice) {
    case GracklePropertyEnum::cooling_time:
      return {&local_calculate_cooling_time,     "cooling_time",     true};
    case GracklePropertyEnum::dust_temperature:
      return {&local_calculate_dust_temperature, "dust_temperature", true};
    case GracklePropertyEnum::gamma:
      return {&local_calculate_gamma,            "gamma",            false};
    case GracklePropertyEnum::pressure:
      return {&local_calculate_pressure,         "pressure",         false};
    case GracklePropertyEnum::temperature:
      return {&local_calculate_temperature,      "temperature",      false};
    default:
      ERROR("GrackleFacade::compute_local_property_", "unknown func_choice");
  }
}
#endif /* CONFIG_USE_GRACKLE */

} // anonymous namespace

//----------------------------------------------------------------------------

void GrackleFacade::compute_property
(const EnzoFieldAdaptor& fadaptor, GracklePropertyEnum func_choice,
 enzo_float* values, int stale_depth, grackle_field_data* grackle_fields)
  const noexcept
{
#ifndef CONFIG_USE_GRACKLE
  ERROR("GrackleFacade::compute_local_property_", "grackle isn't being used");
#else

  GracklePropFnInfo fn_info = get_fn_info_(func_choice);

  // setup code_units struct
  code_units my_units;
  double compute_time =
    (enzo::cosmology() != nullptr) ? fadaptor.compute_time() : -1.0;
  setup_grackle_u_(compute_time, radiation_redshift_, &my_units);

  grackle_field_data cur_grackle_fields;

  // if grackle fields are not provided, define them
  bool delete_grackle_fields = false;
  if (!grackle_fields){
    // the cell width is not used for computing for local properties. Thus, we
    // don't require it (for cases where fadaptor wraps EnzoEFltArrayMap)
    bool omit_cell_width = true;

    grackle_fields  = &cur_grackle_fields;
    this->setup_grackle_fields(fadaptor, grackle_fields,
                               stale_depth, omit_cell_width);
    delete_grackle_fields = true;
  }

  for (int i = 1; i <= grackle_fields->grid_rank; i++){
    int ax_start = grackle_fields->grid_start[i-1];
    int ax_end = grackle_fields->grid_end[i-1];
    int ax_dim = grackle_fields->grid_dimension[i-1];

    // previously Grackle's local_calculate_pressure, local_calculate_gamma,
    // & local_calculate_temperature functions ignored grid_start & grid_end.
    // This should no longer be a problem (it was fixed in PR #106 and we
    // require a more recent version of grackle), but currently this is
    // untested
    ASSERT("GrackleFacade::compute_local_property_",
           ("this method is untested when grackle_fields->grid_start & "
            "grackle_fields->grid_end don't include all data."),
           (ax_start == 0) & ((ax_end+1) == ax_dim));
  }

  // because this function is const-qualified, my_chemistry_.get_ptr() and
  // grackle_rates_ are currently pointers to constants
  // we need to drop the `const` to be able to pass to these values to func
  // (this is okay because func will not actually modify the value).
  chemistry_data * chemistry_data_ptr
    = const_cast<chemistry_data *>(my_chemistry_.get_ptr());
  chemistry_data_storage * grackle_rates_ptr
    = const_cast<chemistry_data_storage *>(grackle_rates_.get());

  if ((*(fn_info.ptr))(chemistry_data_ptr, grackle_rates_ptr,
                       &my_units, grackle_fields, values) == ENZO_FAIL){
    ERROR1("GrackleFacade::compute_local_property_()",
	   "Error during call to Grackle's local_calculate_%s routine",
           fn_info.name_suffix);
  }
  if (delete_grackle_fields){
    this->delete_grackle_fields(grackle_fields);
  }
#endif /* CONFIG_USE_GRACKLE */
}
