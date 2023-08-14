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
  ///
  /// @todo consider removing this constructor (it's never used)
  EnzoComputeTemperature
  (double density_floor,
   double temperature_floor,
   double mol_weight,
   bool comoving_coordinates);

  /// Create a new EnzoComputeTemperature object
  EnzoComputeTemperature(const EnzoPhysicsFluidProps* fluid_props,
                         bool comoving_coordinates)
    : EnzoComputeTemperature(0.0, 0.0, 0.0, comoving_coordinates)
  {
    ASSERT("EnzoComputeTemperature::EnzoComputeTemperature",
           "does not accept a nullptr", fluid_props != nullptr);
    const EnzoFluidFloorConfig& floor_conf = fluid_props->fluid_floor_config();

    if (!fluid_props->eos_variant().holds_alternative<EnzoEOSIdeal>()){
      ERROR("EnzoComputeTemperature::EnzoComputeTemperature",
            "Currently unsure how to compute the temperature for the current "
            "type of EOS.");
    }

    density_floor_ = floor_conf.has_density_floor() ?
      (double)floor_conf.density() : 0.0;
    temperature_floor_ = floor_conf.has_temperature_floor() ?
      (double)floor_conf.temperature() : 0.0;
    mol_weight_ = (double)fluid_props->mol_weight();
  }

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
  ///
  /// This recomputes the "pressure" field and overwrites the values
  virtual void compute( Block * block) throw();

  virtual void compute( Block * block, enzo_float * t) throw();

  // name of derived field that this function calculates
  std::string name () throw() {
    return "temperature";
  }

  void compute_(Block * block,
    enzo_float * t,
    bool recompute_presure = true,
    grackle_field_data * grackle_fields = nullptr
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
