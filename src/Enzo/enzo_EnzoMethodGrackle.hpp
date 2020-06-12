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
    p | time_grackle_data_initialized_;

    if (p.isUnpacking()) {
      // the following recomputes grackle_rates_. This avoids having to write
      // pup methods for all of Grackle's internal data structures
      initialize_grackle_chemistry_data(time_grackle_data_initialized_);
    }

  #endif /* CONFIG_USE_GRACKLE */

  }

  /// Apply the method to advance a block one timestep
  virtual void compute( Block * block) throw();

  virtual std::string name () throw ()
  { return "grackle"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

#ifdef CONFIG_USE_GRACKLE

  void initialize_grackle_chemistry_data(double current_time);

  static void setup_grackle_units(EnzoBlock * enzo_block,
                                  code_units * grackle_units,
                                  int i_hist = 0 ) throw();

  static void setup_grackle_fields(EnzoBlock * enzo_block,
                                   grackle_field_data * grackle_fields,
                                   int i_hist = 0 ) throw();

  static void update_grackle_density_fields(EnzoBlock * enzo_block,
                                   grackle_field_data * grackle_fields) throw();

  static void delete_grackle_fields(grackle_field_data * grackle_fields_) throw() {
      grackle_fields_->density         = NULL;
      grackle_fields_->internal_energy = NULL;
      grackle_fields_->x_velocity      = NULL;
      grackle_fields_->y_velocity      = NULL;
      grackle_fields_->z_velocity      = NULL;
      grackle_fields_->HI_density      = NULL;
      grackle_fields_->HII_density     = NULL;
      grackle_fields_->HeI_density     = NULL;
      grackle_fields_->HeII_density    = NULL;
      grackle_fields_->HeIII_density   = NULL;
      grackle_fields_->e_density       = NULL;
      grackle_fields_->HM_density      = NULL;
      grackle_fields_->H2I_density     = NULL;
      grackle_fields_->H2II_density    = NULL;
      grackle_fields_->DI_density      = NULL;
      grackle_fields_->DII_density     = NULL;
      grackle_fields_->HDI_density     = NULL;
      grackle_fields_->metal_density   = NULL;
      grackle_fields_->volumetric_heating_rate = NULL;
      grackle_fields_->specific_heating_rate   = NULL;

      delete [] grackle_fields_->grid_dimension; grackle_fields_->grid_dimension = NULL;
      delete [] grackle_fields_->grid_start;     grackle_fields_->grid_start      = NULL;
      delete [] grackle_fields_->grid_end;       grackle_fields_->grid_end        = NULL;

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
