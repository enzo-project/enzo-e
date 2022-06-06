// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSIdeal.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of the equation of state for an ideal
/// adiabatic gas

#ifndef ENZO_ENZO_EOS_IDEAL_HPP
#define ENZO_ENZO_EOS_IDEAL_HPP

class EnzoEOSIdeal : public EnzoEquationOfState
{

  /// @class    EnzoEOSIdeal
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates equation of state for ideal gas
  
public: // interface
  
  /// Create a new EnzoEOSIdeal object
  EnzoEOSIdeal(double gamma, double density_floor, double pressure_floor,
	       bool dual_energy_formalism,
	       double dual_energy_formalism_eta) throw()
    : EnzoEquationOfState(),
      gamma_(gamma),
      density_floor_(density_floor),
      pressure_floor_(pressure_floor),
      dual_energy_formalism_(dual_energy_formalism),
      dual_energy_formalism_eta_(dual_energy_formalism_eta)
  { }

  /// Delete EnzoEOSIdeal object
  ~EnzoEOSIdeal()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoEOSIdeal);

  /// CHARM++ migration constructor for PUP::able
  EnzoEOSIdeal (CkMigrateMessage *m)
    : EnzoEquationOfState(m),
      gamma_(0.),
      density_floor_(0.),
      pressure_floor_(0.),
      dual_energy_formalism_(false),
      dual_energy_formalism_eta_(0.)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void primitive_from_integration
  (const EnzoEFltArrayMap &integration_map, EnzoEFltArrayMap &primitive_map,
   const int stale_depth, const str_vec_t &passive_list) const;

  void pressure_from_integration
  (const EnzoEFltArrayMap &integration_map,
   const CelloArray<enzo_float, 3> &pressure,
   const int stale_depth) const;

  inline enzo_float get_density_floor() const { return density_floor_; }

  enzo_float get_pressure_floor() const { return pressure_floor_; }

  void apply_floor_to_energy_and_sync(EnzoEFltArrayMap &integration_map,
                                      const int stale_depth) const;

  bool is_barotropic() const { return false; }

  enzo_float get_gamma() const { return gamma_;}

  enzo_float get_isothermal_sound_speed() const { return 0;}

  // In the future, this won't be hardcoded to false
  bool uses_dual_energy_formalism() const { return dual_energy_formalism_; };

protected: // attributes
  enzo_float gamma_; // adiabatic index
  enzo_float density_floor_;
  enzo_float pressure_floor_;
  bool dual_energy_formalism_;
  double dual_energy_formalism_eta_;
};

#endif /* ENZO_ENZO_EOS_IDEAL_HPP */
