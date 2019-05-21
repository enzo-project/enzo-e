// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputeTemperature.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
///           Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-05-07
/// @brief    [\ref Enzo] Implementation of Enzo's ComputeTemperature functions

#ifndef ENZO_ENZO_COMPUTE_TEMPERATURE_HPP
#define ENZO_ENZO_COMPUTE_TEMPERATURE_HPP

class EnzoComputeTemperature : public Compute {

  /// @class    EnzoComputeTemperature
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ComputeTemperature functions

public: // interface

  /// Create a new EnzoComputeTemperature object
  EnzoComputeTemperature
  (double density_floor,
   double temperature_floor,
   double mol_weight,
   bool comoving_coordinates);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoComputeTemperature);

  /// Charm++ PUP::able migration constructor
  EnzoComputeTemperature (CkMigrateMessage *m)
    : Compute(m),
      density_floor_(0.0),
      temperature_floor_(0.0),
      mol_weight_(0.0),
      comoving_coordinates_(false)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Perform the computation on the block
  virtual void compute( Block * block) throw();

  virtual void compute( Block * block, enzo_float * t) throw();

  // name of derived field that this function calculates
  std::string name () throw() {
    return "temperature";
  }

  void compute_(Block * block,
    enzo_float * t,
    bool recompute_presure = true
#ifdef CONFIG_USE_GRACKLE
 , code_units * grackle_units = NULL,
   grackle_field_data * grackle_fields = NULL
#endif
 );

private: // functions


private: // attributes

  // minimum density
  double density_floor_;

  // minimum temperature: default 1.0
  double temperature_floor_;

  // mol weight: default 0.6
  double mol_weight_;

  bool comoving_coordinates_;

};

#endif /* ENZO_ENZO_COMPUTE_TEMPERATURE_HPP */
