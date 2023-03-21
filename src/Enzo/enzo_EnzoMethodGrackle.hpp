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

class EnzoMethodGrackle : public Method {

  /// @class    EnzoMethodGrackle
  /// @ingroup  Enzo
  ///
  /// This class interfaces the Grackle primordial chemistry / cooling
  /// library with Cello
  ///
  /// In more detail, this class implements higher-level functionallity,
  /// (e.g. related to the intricacies of implementing a ``Method`` subclass),
  /// while the ``GrackleFacade`` attribute is responsible for low-level
  /// Grackle details (e.g. managing the lifetime of Grackle's configuration
  /// structs).

public: // interface

  /// Create a new EnzoMethodGrackle object
  EnzoMethodGrackle(GrackleChemistryData my_chemistry,
                    bool use_cooling_timestep,
                    const double radiation_redshift,
                    const double physics_cosmology_initial_redshift,
                    const double time);

  // Destructor
  virtual ~EnzoMethodGrackle() throw() { }

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGrackle);

  /// Charm++ PUP::able migration constructor
  EnzoMethodGrackle (CkMigrateMessage *m)
    : Method (m),
      grackle_facade_(m),
      use_cooling_timestep_(false)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    // NOTE: change this function whenever attributes change

    TRACEPUP;

    Method::pup(p);
    p | grackle_facade_;
    p | use_cooling_timestep_;
  }

  /// Apply the method to advance a block one timestep
  virtual void compute( Block * block) throw();

  virtual std::string name () throw ()
  { return "grackle"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) throw();

  /// returns the stored instance of GrackleChemistryData, if the simulation is
  /// configured to actually use grackle
  ///
  /// In practice, if this were to ever return a nullptr, then this
  /// EnzoMethodGrackle instance probably shouldn't exist (the program will
  /// probably abort somewhere along the lines)
  const GrackleChemistryData* try_get_chemistry() const throw() {
    return grackle_facade_.try_get_chemistry();
  }

  void define_required_grackle_fields();

  void setup_grackle_fields(const EnzoFieldAdaptor& fadaptor,
                            grackle_field_data * grackle_fields,
                            int stale_depth = 0,
                            bool omit_cell_width = false) const throw()
  {
    grackle_facade_.setup_grackle_fields(fadaptor, grackle_fields,
                                         stale_depth, omit_cell_width);
  }

  void setup_grackle_fields(Block * block,
                            grackle_field_data * grackle_fields,
                            int i_hist = 0 ) const throw()
  {
    grackle_facade_.setup_grackle_fields(block, grackle_fields, i_hist);
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
    grackle_facade_.delete_grackle_fields(grackle_fields);
  }


  void enforce_metallicity_floor(Block * block) throw();

  void calculate_cooling_time(const EnzoFieldAdaptor& fadaptor,
                              enzo_float* ct, int stale_depth = 0,
			      grackle_field_data* grackle_fields = nullptr)
    const throw()
  {
    grackle_facade_.compute_property(fadaptor,
                                     GracklePropertyEnum::cooling_time,
                                     ct, stale_depth, grackle_fields);
  }

  void calculate_pressure(const EnzoFieldAdaptor& fadaptor,
                          enzo_float* pressure, int stale_depth = 0,
			  grackle_field_data* grackle_fields = nullptr)
    const throw()
  {
    grackle_facade_.compute_property(fadaptor, GracklePropertyEnum::pressure,
                                     pressure, stale_depth, grackle_fields);
  }

  void calculate_temperature(const EnzoFieldAdaptor& fadaptor,
                             enzo_float* temperature, int stale_depth = 0,
			     grackle_field_data* grackle_fields = nullptr)
    const throw()
  {
    grackle_facade_.compute_property(fadaptor, GracklePropertyEnum::temperature,
                                     temperature, stale_depth, grackle_fields);
  }

protected: // methods

  void compute_( Block * block) throw();

  void ResetEnergies ( Block * block) throw();

protected: // attributes
  /// the GrackleFacade instance provides an interface to all operations in the
  /// Grackle library and stores the current configuration. You can assume that
  /// this is always correctly initialized
  GrackleFacade grackle_facade_;
  bool use_cooling_timestep_;
};

#endif /* ENZO_ENZO_METHOD_GRACKLE_HPP */
