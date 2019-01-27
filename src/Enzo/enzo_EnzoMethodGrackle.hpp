// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoMethodGrackle.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     Thu May 15 14:32:28 EDT 2014
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
                    const float physics_cosmology_initial_redshift,
                    const float time);

  // Destructor
  virtual ~EnzoMethodGrackle() throw() {}

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodGrackle);

  /// Charm++ PUP::able migration constructor
  EnzoMethodGrackle (CkMigrateMessage *m)
    : Method (m)
#ifdef CONFIG_USE_GRACKLE
      , grackle_units_()
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

  #endif /* CONFIG_USE_GRACKLE */

  }

  /// Apply the method to advance a block one timestep
  virtual void compute( Block * block) throw();

  virtual std::string name () throw ()
  { return "grackle"; }

  /// Compute maximum timestep for this method
  virtual double timestep ( Block * block) const throw();

#ifdef CONFIG_USE_GRACKLE
  static void setup_grackle_units(EnzoBlock * enzo_block,
                                  code_units * grackle_units) throw();

  static void setup_grackle_fields(EnzoBlock * enzo_block,
                                   grackle_field_data * grackle_fields,
                                   bool initialize_fields = false) throw();

  static void update_grackle_density_fields(EnzoBlock * enzo_block,
                                   grackle_field_data * grackle_fields) throw();
#endif

protected: // methods

#ifdef CONFIG_USE_GRACKLE
  void compute_( EnzoBlock * enzo_block) throw();

  void ResetEnergies ( EnzoBlock * enzo_block) throw();

// protected: // attributes

  code_units grackle_units_;
#endif /* ENZO_ENZO_METHOD_GRACKLE_HPP */


};

#endif /* ENZO_ENZO_METHOD_GRACKLE_HPP */
