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
      temperature_(0.0)
  {}

  EnzoFluidFloorConfig(enzo_float density_floor,
                       enzo_float pressure_floor,
                       enzo_float temperature_floor)
    : density_(density_floor),
      pressure_(pressure_floor),
      temperature_(temperature_floor)
  {}

  bool has_density_floor() const noexcept { return density_ > 0; }
  bool has_pressure_floor() const noexcept { return pressure_ > 0; }
  bool has_temperature_floor() const noexcept {return temperature_ > 0; }

  /// query the density floor. The program aborts if it is not set.
  enzo_float density() const noexcept
  {
    ASSERT("EnzoFluidFloorConfig::density",
           "this instance does not hold a density floor",
           has_density_floor());
    return density_;
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

  void pup(PUP::er &p) {
    p|density_;
    p|pressure_;
    p|temperature_;
  }

private: // attributes

  enzo_float density_;
  enzo_float pressure_;
  enzo_float temperature_;
};

#endif /* ENZO_ENZO_FLUID_FLOORS_HPP */
