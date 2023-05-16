// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoDualEnergyConfig.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2022-04-05
/// @brief    [\ref Enzo] Declaration of the EnzoDualEnergyConfig class

#ifndef ENZO_ENZO_DUAL_ENERGY_CONFIG_HPP
#define ENZO_ENZO_DUAL_ENERGY_CONFIG_HPP

class EnzoDualEnergyConfig {

  /// @class    EnzoDualEnergyConfig
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the configuration of the dual-energy
  ///    formalism
  ///
  /// This can represent 3 states:
  /// 1. "disabled" - the dual energy is not in use
  /// 2. "modern" - holds a single parameter, eta. For more details see the
  ///    description on the website of the dual-energy formalism that is used
  ///    by the VL+CT integrator
  /// 2. "bryan95" - holds the parameters, eta1 and eta2, to be used according
  ///    to the procedure described in bryan+ (1995)
  ///    https://ui.adsabs.harvard.edu/abs/1995CoPhC..89..149B/abstract
  ///
  /// The data stored in the class should be considered immutable.

public:

  /// creates an instance that represents a disabled dual energy formalism
  EnzoDualEnergyConfig()
    : EnzoDualEnergyConfig(-1, -1)
  { }

  /// creates an instance that represents a disabled dual energy formalism
  static EnzoDualEnergyConfig build_disabled() noexcept
  { return EnzoDualEnergyConfig(); }

  /// creates an instance that holds parameters for the modern formulation of
  /// the dual energy formalism
  static EnzoDualEnergyConfig build_modern_formulation(enzo_float eta) noexcept
  {
    ASSERT("EnzoDualEnergyConfig::build_modern_formulation",
           "eta must be non-negative", eta >= 0);
    return {eta, -1};
  }

  /// creates an instance that holds parameters for the bryan95 formulation of
  /// the dual energy formalism
  static EnzoDualEnergyConfig build_bryan95_formulation
  (enzo_float eta1, enzo_float eta2) noexcept
  {
    ASSERT("EnzoDualEnergyConfig::build_bryan95_formulation",
           "eta1 and eta2 must be non-negative", (eta1 >= 0) & (eta2 >= 0));
    return {eta1, eta2};
  }

  /// queries if any formulation of the dual energy formalism is enabled
  bool any_enabled() const noexcept { return !is_disabled(); }

  /// queries if the dual energy formalism is disabled
  bool is_disabled() const noexcept { return primary_eta_ < 0; }

  /// queries whether this instance holds configuration parameters for the
  /// modern formulation of the dual energy formalism (and can be used to
  /// access the values of the parameters)
  ///
  /// @param[out] eta When this instance holds values for the modern
  ///    formulation of the dual-energy formalism, the corresponding eta value
  ///    parameter is written to the location specified by this pointer (unless
  ///    the parameter is a nullptr).
  bool modern_formulation(enzo_float* eta = nullptr) const noexcept
  {
    if ((primary_eta_ >= 0) & (other_eta_ < 0)){
      if (eta != nullptr) { *eta = primary_eta_; }
      return true;
    } else {
      return false;
    }
  }

  /// queries whether this instance holds configuration parameters for the
  /// bryan95 formulation of the dual energy formalism (and can be used to
  /// access the values of the parameters)
  ///
  /// @param[out] eta1, eta2 When this instance holds values for the bryan95
  ///    formulation of the dual-energy formalism, the corresponding parameter
  ///    values are written to locations specified by these pointers (unless
  ///    the parameters are nullptrs)
  bool bryan95_formulation(enzo_float* eta1 = nullptr,
                           enzo_float* eta2 = nullptr) const noexcept
  {
    if ((primary_eta_ >= 0) & (other_eta_ >= 0)){
      if (eta1 != nullptr) { *eta1 = primary_eta_; }
      if (eta2 != nullptr) { *eta2 = other_eta_; }
      return true;
    } else {
      return false;
    }
  }

  void pup(PUP::er &p) {
    p|primary_eta_;
    p|other_eta_;
  }
  
private:

  /// constructs a new EnzoDualEnergyConfig instance. Users don't directly use
  /// this method, instead they should should use one of the static "build_*"
  /// methods instead.
  EnzoDualEnergyConfig(enzo_float primary_eta, enzo_float other_eta)
    : primary_eta_(primary_eta),
      other_eta_(other_eta)
  {
    if ((primary_eta < 0) & (other_eta >= 0)){
      ERROR("EnzoDualEnergyConfig::EnzoDualEnergyConfig", "invalid arguments");
    }
  }

protected: // attributes

  enzo_float primary_eta_;
  enzo_float other_eta_;
};

#endif /* ENZO_ENZO_DUAL_ENERGY_CONFIG_HPP */
