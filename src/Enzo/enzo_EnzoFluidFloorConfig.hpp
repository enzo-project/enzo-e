// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoFluidFloorConfig.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-05-01
/// @brief    [\ref Enzo] Declaration of the EnzoFluidFloorConfig class

#ifndef ENZO_ENZO_FLUID_FLOOR_CONFIG_HPP
#define ENZO_ENZO_FLUID_FLOOR_CONFIG_HPP

class EnzoFluidFloorConfig {

  /// @class    EnzoFluidFloorConfig
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Stores the floors for various fluid quantities
  ///
  /// The data stored in the class should be considered immutable.
  ///
  /// As of now a valid floor must be positive.
  ///
  /// Going forward, we may want to revisit this class's implementation.

public: // interface

  EnzoFluidFloorConfig()
    : density_(0.0),
      pressure_(0.0),
      temperature_(0.0),
      metal_mass_frac_(0.0)
  {}

  EnzoFluidFloorConfig(double density_floor,
                       double pressure_floor,
                       double temperature_floor,
                       double metal_mass_frac)
    : density_(density_floor),
      pressure_(pressure_floor),
      temperature_(temperature_floor),
      metal_mass_frac_(metal_mass_frac)
  {}

  bool has_density_floor() const noexcept { return density_ > 0; }
  bool has_pressure_floor() const noexcept { return pressure_ > 0; }
  bool has_temperature_floor() const noexcept {return temperature_ > 0; }
  bool has_metal_mass_frac_floor() const noexcept
  { return metal_mass_frac_ > 0; }

  /// query the density floor. The value is always returned in double precision
  ///
  /// This is only present as a short-term solution. We should remove this
  /// method in the near future.
  double density_dbl_prec() const noexcept
  {
    ASSERT("EnzoFluidFloorConfig::density_prec",
           "this instance does not hold a density floor",
           has_density_floor());
    return density_;
  }

  /// query the density floor. The program aborts if it is not set.
  enzo_float density() const noexcept
  {
    ASSERT("EnzoFluidFloorConfig::density",
           "this instance does not hold a density floor",
           has_density_floor());
    return static_cast<enzo_float>(density_);
  }

  /// query the pressure floor. The program aborts if it is not set.
  enzo_float pressure() const noexcept
  {
    ASSERT("EnzoFluidFloorConfig::pressure",
           "this instance does not hold a pressure floor",
           has_pressure_floor());
    return pressure_;
  }

  /// query the temperature floor. The program aborts if it is not set.
  enzo_float temperature() const noexcept
  {
    ASSERT("EnzoFluidFloorConfig::temperature",
           "this instance does not hold a temperature floor",
           has_temperature_floor());
    return temperature_;
  }

  /// query the metal mass fraction floor. The program aborts if it is not set.
  enzo_float metal_mass_frac() const noexcept
  {
    ASSERT("EnzoFluidFloorConfig::metal_mass_frac",
           "this instance does not hold a floor for the metal mass fraction",
           has_metal_mass_frac_floor());
    return metal_mass_frac_;
  }

  void pup(PUP::er &p) {
    p|density_;
    p|pressure_;
    p|temperature_;
    p|metal_mass_frac_;
  }

private: // attributes

  double density_;
  enzo_float pressure_;
  enzo_float temperature_;
  enzo_float metal_mass_frac_;
};

#endif /* ENZO_ENZO_FLUID_FLOORS_HPP */
