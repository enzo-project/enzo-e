// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSIsothermal.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-04-01
/// @brief    [\ref Enzo] Declaration of the EnzoEOSIsothermal struct

#ifndef ENZO_ENZO_EOS_ISOTHERMAL_HPP
#define ENZO_ENZO_EOS_ISOTHERMAL_HPP

struct EnzoEOSIsothermal {
  /// @class    EnzoEOSIsothermal
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the equation of state for an
  ///           isothermal gas.
  ///
  /// At the moment, this mostly just acts as a placeholder that serves 2
  /// purposes:
  ///    1. To act as a dummy alternative to the EnzoEOSIdeal
  ///    2. To act as the EOS object that should be active while using the PPML
  ///       solver

public: // attributes:

  // currently there are no attributes (in the future, we might store stuff
  // like isothermal sound-speed)

public: // public interface common to all EOS types
  constexpr static const char* name() noexcept { return "isothermal"; }

  constexpr static bool is_barotropic() noexcept { return true; }

  /// create a string for debugging purposes that represents the EOS's value
  std::string debug_string() const noexcept { return "EnzoEOSIsothermal{}"; }

};

/// function responsible for PUPing
///
/// We have explicitly chosen not to making this a method of the EOS class
/// to try to avoid issues with copying instances to a GPU
inline void pup(PUP::er &p, EnzoEOSIsothermal& eos) noexcept {
  // make sure to keep this synchronized with the attributes EnzoEOSIdeal

  // currently, there are no attributes to pup
}

#endif /* ENZO_ENZO_EOS_ISOTHERMAL_HPP */
