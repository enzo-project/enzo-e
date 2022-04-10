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
  /// static method to compute pressure
  static void compute_pressure(const EnzoFieldAdaptor& fadaptor,
                               const CelloArray<enzo_float, 3>& p,
                               bool mhd, bool dual_energy, double gamma,
                               int stale_depth = 0
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
