// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoComputePressure.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
///           Andrew Emerick (aemerick11@gmail.com)
/// @date     2019-05-07
/// @brief    [\ref Enzo] Implementation of Enzo's ComputePressure functions

#ifndef ENZO_ENZO_COMPUTE_PRESSURE_HPP
#define ENZO_ENZO_COMPUTE_PRESSURE_HPP

class EnzoComputePressure : public Compute {

  /// @class    EnzoComputePressure
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Enzo's ComputePressure functions

public: // interface

  /// Create a new EnzoComputePressure object
  EnzoComputePressure(double gamma,
		      bool comoving_coordinates);

  /// Charm++ PUP::able declarations
  PUPable_decl(EnzoComputePressure);

  /// Charm++ PUP::able migration constructor
  EnzoComputePressure (CkMigrateMessage *m)
    : Compute(m),
      gamma_(0.0),
      comoving_coordinates_(false)
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  // name of derived field that this function calculates
  std::string name () throw() {
    return "pressure";
  }

  /// Perform the computation on the block and store the results in the
  /// "pressure" field
  void compute( Block * block) throw();

  /// Perform the computation on the block and store the result in the provided
  /// array.
  ///
  /// @param[in]  block Where input arrays are loaded from
  /// @param[out] p Array that should hold the result.
  /// @param[in] stale_depth The number of the outermost zones to exclude from
  ///     the calculation
  ///
  /// This function will rely on the Grackle-supplied function for computing
  /// the pressure if the simulation is configured to use `EnzoMethodGrackle`
  /// (in this case the `gamma` value passed to the constructor is ignored).
  /// Otherwise, it falls back to an alternate version of the calculation.
  void compute( Block * block,
                enzo_float* p,
                int stale_depth = 0) throw();

  /// Perform the calculation with the provided `grackle_units` or
  /// `grackle_fields` object.
  ///
  /// @note
  /// This is primarily used by `EnzoInitialGrackleTest`
  void compute_(Block * block,
                enzo_float * p,
                int stale_depth = 0
#ifdef CONFIG_USE_GRACKLE
                , code_units * grackle_units = nullptr,
                grackle_field_data * grackle_fields = nullptr
#endif
                );

  /// static method to compute thermal pressure
  ///
  /// @param[in]  fadaptor Contains arrays of quantities used to compute the
  ///     thermal pressure.
  /// @param[out] p Array where the thermal pressure is to be stored
  /// @param[in]  mhd Whether the simulation includes magnetic fields
  /// @param[in]  dual_energy Whether the simulation uses the dual-energy
  ///     formalism
  /// @param[in]  gamma Specifies the adiabatic index. This is entirely ignored
  ///     if the Grackle-supplied function for computing pressure is used.
  /// @param[in]  stale_depth The number of the outermost zones to exclude from
  ///     the calculation.
  /// @param[in]  ignore_grackle indicates whether to ignore the Grackle
  ///     routine for computing pressure. This exists to force the usage of the
  ///     `gamma` parameter (and can be used to ignore the effects of molecular
  ///     hydrogen on the adiabatic index)
  ///
  /// By default, this function will rely on the Grackle-supplied function for
  /// computing the pressure if the simulation is configured to use
  /// `EnzoMethodGrackle`. The `ignore_grackle` parameter only has meaning in
  /// this scenario (if the simulation does not use `EnzoMethodGrackle`, then
  /// the value of `ignore_grackle` is meaningless).
  static void compute_pressure(const EnzoFieldAdaptor& fadaptor,
                               const CelloArray<enzo_float, 3>& p,
                               bool mhd, bool dual_energy, double gamma,
                               int stale_depth = 0,
                               bool ignore_grackle = false
#ifdef CONFIG_USE_GRACKLE
                              , code_units * grackle_units = nullptr,
                               grackle_field_data * grackle_fields = nullptr
#endif
                               ) throw();

protected: // attributes

  double gamma_;
  bool comoving_coordinates_;

};

#endif /* ENZO_ENZO_COMPUTE_PRESSURE_HPP */
