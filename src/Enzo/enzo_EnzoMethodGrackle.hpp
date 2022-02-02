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
  EnzoMethodGrackle(
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
      , grackle_units_()
      , grackle_rates_()
      , time_grackle_data_initialized_(ENZO_FLOAT_UNDEFINED)
#endif
    {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {

    // NOTE: change this function whenever attributes change

  #ifdef CONFIG_USE_GRACKLE

    TRACEPUP;

    Method::pup(p);

    p | grackle_units_;

    double last_init_time = time_grackle_data_initialized_;
    p | last_init_time;
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

#ifdef CONFIG_USE_GRACKLE

  void define_required_grackle_fields();

  void initialize_grackle_chemistry_data(double current_time,
                                         bool preinitialized_units = false);

  static void setup_grackle_units (double current_time,
                                   code_units * grackle_units) throw();

  static void setup_grackle_units(EnzoBlock * enzo_block,
                                  code_units * grackle_units,
                                  int i_hist = 0 ) throw()
  {
    double compute_time;
    if (i_hist == 0) {
      compute_time = enzo_block->time();
    } else {
      compute_time = enzo_block->data()->field().history_time(i_hist);
    }
    setup_grackle_units(compute_time, grackle_units);
  }

  static void setup_grackle_fields(EnzoBlock * enzo_block,
                                   grackle_field_data * grackle_fields,
                                   int i_hist = 0 ) throw();

  static void update_grackle_density_fields(EnzoBlock * enzo_block,
                                   grackle_field_data * grackle_fields) throw();

  static void delete_grackle_fields(grackle_field_data * grackle_fields) throw() {
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

  void calculate_cooling_time(Block * block, enzo_float* ct,
			      code_units* grackle_units = NULL,
			      grackle_field_data* grackle_fields = NULL,
			      int i_hist = 0) const throw()
  {
    compute_local_property_(block, ct, grackle_units, grackle_fields, i_hist,
			    &local_calculate_cooling_time,
			    "local_calculate_cooling_time");
  }

  void calculate_pressure(Block * block, enzo_float* pressure,
			  code_units* grackle_units = NULL,
			  grackle_field_data* grackle_fields = NULL,
			  int i_hist = 0) const throw()
  {
    compute_local_property_(block, pressure, grackle_units, grackle_fields,
			    i_hist, &local_calculate_pressure,
			    "local_calculate_pressure");
  }

  void calculate_temperature(Block * block, enzo_float* temperature,
			     code_units* grackle_units = NULL,
			     grackle_field_data* grackle_fields = NULL,
			     int i_hist = 0) const throw()
  {
    compute_local_property_(block, temperature, grackle_units, grackle_fields,
			    i_hist, &local_calculate_temperature,
			    "local_calculate_temperature");
  }

#endif

protected: // methods

  void deallocate_grackle_rates_() throw();

#ifdef CONFIG_USE_GRACKLE

  // when grackle_units is NULL, new values are temporarily allocated
  void compute_local_property_(Block * block, enzo_float* values,
			       code_units* grackle_units,
			       grackle_field_data* grackle_fields, int i_hist,
			       grackle_local_property_func func,
			       std::string func_name) const throw();

#endif

protected: // methods

#ifdef CONFIG_USE_GRACKLE
  void compute_( EnzoBlock * enzo_block) throw();

  void ResetEnergies ( EnzoBlock * enzo_block) throw();

// protected: // attributes

  code_units grackle_units_;
  chemistry_data_storage grackle_rates_;
  double time_grackle_data_initialized_;
#endif

};

#endif /* ENZO_ENZO_METHOD_GRACKLE_HPP */
