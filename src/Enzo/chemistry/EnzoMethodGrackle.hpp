// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
///           Andrew Emerick (aemerick11@gmail.com)
/// @date     Tues May 7
/// @brief    [\ref Enzo] Declaration of EnzoMethodGrackle class
///
/// This class interfaces the Grackle primordial chemistry / cooling
/// library with Cello

#ifndef ENZO_ENZO_METHOD_GRACKLE_HPP
#define ENZO_ENZO_METHOD_GRACKLE_HPP


#ifdef CONFIG_USE_GRACKLE
// PUP operator for Grackle's code_units
inline void operator|(PUP::er &p, code_units &c){
  // Make sure to change this if code_units ever changes
  // all are just single values (int or float)

  p | c.comoving_coordinates;
  p | c.density_units;
  p | c.length_units;
  p | c.time_units;
  p | c.velocity_units;
  p | c.a_units;
  p | c.a_value;

}

typedef int (*grackle_local_property_func)(chemistry_data*,
					   chemistry_data_storage *,
					   code_units*, grackle_field_data*,
					   enzo_float*);
#endif



class EnzoMethodGrackle : public Method {

  /// @class    EnzoMethodGrackle
  /// @ingroup  Enzo
  ///
  /// This class interfaces the Grackle primordial chemistry / cooling
  /// library with Cello

public: // interface

  /// Create a new EnzoMethodGrackle object
  EnzoMethodGrackle(GrackleChemistryData my_chemistry,
                    bool use_cooling_timestep,
                    const double radiation_redshift,
                    const double physics_cosmology_initial_redshift,
                    const double time);

  // Destructor
  virtual ~EnzoMethodGrackle() throw() { deallocate_grackle_rates_(); }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGrackle);

  /// Charm++ PUP::able migration constructor
  EnzoMethodGrackle (CkMigrateMessage *m)
    : Method (m)
#ifdef CONFIG_USE_GRACKLE
      , my_chemistry_()
      , grackle_units_()
      , grackle_rates_()
      , time_grackle_data_initialized_(ENZO_FLOAT_UNDEFINED)
      , radiation_redshift_(-1.0)
      , use_cooling_timestep_(false)
#endif
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    // NOTE: change this function whenever attributes change

  #ifdef CONFIG_USE_GRACKLE

    TRACEPUP;

    Method::pup(p);
    p | my_chemistry_;
    p | grackle_units_;
    double last_init_time = time_grackle_data_initialized_;
    p | last_init_time;
    p | radiation_redshift_;
    p | use_cooling_timestep_;
    if (p.isUnpacking()) {
      ASSERT("EnzoMethodGrackle::pup",
             "grackle_chemistry_data must have previously been initialized",
             last_init_time!=ENZO_FLOAT_UNDEFINED);
      // the following recomputes grackle_rates_ (and sets the value of
      // time_grackle_data_initialized_ to last_init_time). This is done
      // to avoid writing pup methods for all of Grackle's internal data
      // structures.
      time_grackle_data_initialized_ = ENZO_FLOAT_UNDEFINED;
      initialize_grackle_chemistry_data(last_init_time, true);
    }

  #endif /* CONFIG_USE_GRACKLE */

  }

  /// Apply the method to advance a block one timestep
  virtual void compute( Block * block) throw();

  virtual std::string name () throw ()
  { return "grackle"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();

  /// returns the stored instance of GrackleChemistryData, if the simulation is
  /// configured to actually use grackle
  const GrackleChemistryData* try_get_chemistry() const throw() {
#ifndef CONFIG_USE_GRACKLE
    return nullptr;
#else
    return (my_chemistry_.get<int>("use_grackle") == 1) ?
      &my_chemistry_ : nullptr;
#endif
  }

#ifdef CONFIG_USE_GRACKLE

  void define_required_grackle_fields();

  void initialize_grackle_chemistry_data(double current_time,
                                         bool preinitialized_units = false);

  /// sets up grackle units
  ///
  /// @param[in]  current_time The current time. A negative value can be passed
  ///     if the current value is not known or convenient to get. In that case,
  ///     the program will abort if the current time is needed (i.e. because
  ///     this is a cosmological simulation)
  /// @param[out] grackle_units The object pointed to by this pointer is set up
  ///
  /// @note
  /// This used to be a global method, but there isn't much need for any other
  /// object to use this method. Some of the Compute object still nominally
  /// allow code_units to be passed as an option, but that really isn't used
  /// anywhere (we plan to drop that in the near future).
  void setup_grackle_units (double current_time,
                            code_units * grackle_units) const throw();

  void setup_grackle_units(const EnzoFieldAdaptor& fadaptor,
                           code_units * grackle_units) const throw();

  void setup_grackle_fields(const EnzoFieldAdaptor& fadaptor,
                            grackle_field_data * grackle_fields,
                            int stale_depth = 0,
                            bool omit_cell_width = false) const throw();

  void setup_grackle_fields(Block * block,
                            grackle_field_data * grackle_fields,
                            int i_hist = 0 ) const throw()
  {
    setup_grackle_fields(EnzoFieldAdaptor(block, i_hist),
                         grackle_fields, 0, false);
  }

  /// Assists with problem initialization
  ///
  /// Scales species density fields to be sensible mass fractions of the
  /// initial density field. Problem types that require finer-tuned control
  /// over individual species fields should adapt this function
  /// in their initialization routines.
  void update_grackle_density_fields
  (Block * block, grackle_field_data * grackle_fields = nullptr)
    const throw();

  void delete_grackle_fields(grackle_field_data * grackle_fields) const throw()
  {
      grackle_fields->density         = NULL;
      grackle_fields->internal_energy = NULL;
      grackle_fields->x_velocity      = NULL;
      grackle_fields->y_velocity      = NULL;
      grackle_fields->z_velocity      = NULL;
      grackle_fields->HI_density      = NULL;
      grackle_fields->HII_density     = NULL;
      grackle_fields->HeI_density     = NULL;
      grackle_fields->HeII_density    = NULL;
      grackle_fields->HeIII_density   = NULL;
      grackle_fields->e_density       = NULL;
      grackle_fields->HM_density      = NULL;
      grackle_fields->H2I_density     = NULL;
      grackle_fields->H2II_density    = NULL;
      grackle_fields->DI_density      = NULL;
      grackle_fields->DII_density     = NULL;
      grackle_fields->HDI_density     = NULL;
      grackle_fields->metal_density   = NULL;
      grackle_fields->volumetric_heating_rate = NULL;
      grackle_fields->specific_heating_rate   = NULL;

      delete [] grackle_fields->grid_dimension; grackle_fields->grid_dimension = NULL;
      delete [] grackle_fields->grid_start;     grackle_fields->grid_start      = NULL;
      delete [] grackle_fields->grid_end;       grackle_fields->grid_end        = NULL;

      return;
  }

  void enforce_metallicity_floor(Block * block) throw();

  void calculate_cooling_time(const EnzoFieldAdaptor& fadaptor, enzo_float* ct,
                              int stale_depth = 0,
			      code_units* grackle_units = nullptr,
			      grackle_field_data* grackle_fields = nullptr)
    const throw()
  {
    compute_local_property_(fadaptor, ct, stale_depth, grackle_units,
                            grackle_fields, &local_calculate_cooling_time,
			    "local_calculate_cooling_time");
  }

  void calculate_pressure(const EnzoFieldAdaptor& fadaptor,
                          enzo_float* pressure, int stale_depth = 0,
			  code_units* grackle_units = nullptr,
			  grackle_field_data* grackle_fields = nullptr)
    const throw()
  {
    compute_local_property_(fadaptor, pressure, stale_depth, grackle_units,
                            grackle_fields, &local_calculate_pressure,
			    "local_calculate_pressure");
  }

  void calculate_temperature(const EnzoFieldAdaptor& fadaptor,
                             enzo_float* temperature, int stale_depth = 0,
			     code_units* grackle_units = nullptr,
			     grackle_field_data* grackle_fields = nullptr)
    const throw()
  {
    compute_local_property_(fadaptor, temperature, stale_depth, grackle_units,
                            grackle_fields, &local_calculate_temperature,
			    "local_calculate_temperature");
  }


#endif


protected: // methods

  void deallocate_grackle_rates_() throw();

#ifdef CONFIG_USE_GRACKLE

  // when grackle_units is NULL, new values are temporarily allocated
  void compute_local_property_(const EnzoFieldAdaptor& fadaptor,
                               enzo_float* values, int stale_depth,
			       code_units* grackle_units,
			       grackle_field_data* grackle_fields,
			       grackle_local_property_func func,
			       std::string func_name) const throw();

#endif

protected: // methods

#ifdef CONFIG_USE_GRACKLE
  void compute_( Block * block) throw();

  void ResetEnergies ( Block * block) throw();

#endif

protected: // attributes
#ifdef CONFIG_USE_GRACKLE

  GrackleChemistryData my_chemistry_;
  code_units grackle_units_;
  chemistry_data_storage grackle_rates_;
  double time_grackle_data_initialized_;

  /// the following parameter is only used in non-cosmological simulations. It
  /// specifies the redshift of the UV background
  double radiation_redshift_;

  bool use_cooling_timestep_;

#endif

};

#endif /* ENZO_ENZO_METHOD_GRACKLE_HPP */
