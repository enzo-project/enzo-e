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
  EnzoEOSIdeal(double gamma) throw()
    : EnzoEquationOfState(),
      gamma_(gamma)
  { }

  /// Delete EnzoEOSIdeal object
  ~EnzoEOSIdeal()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoEOSIdeal);

  /// CHARM++ migration constructor for PUP::able
  EnzoEOSIdeal (CkMigrateMessage *m)
    : EnzoEquationOfState(m),
      gamma_(0.)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void primitive_from_integration
  (const EnzoEFltArrayMap &integration_map, EnzoEFltArrayMap &primitive_map,
   const int stale_depth, const str_vec_t &passive_list,
   const bool ignore_grackle) const;

  void pressure_from_integration
  (const EnzoEFltArrayMap &integration_map,
   const CelloArray<enzo_float, 3> &pressure,
   const int stale_depth, const bool ignore_grackle) const;

  void apply_floor_to_energy_and_sync(EnzoEFltArrayMap &integration_map,
                                      const int stale_depth) const;

  bool is_barotropic() const { return false; }

  enzo_float get_gamma() const { return gamma_;}

  enzo_float get_isothermal_sound_speed() const { return 0;}

protected: // attributes
  enzo_float gamma_; // adiabatic index
};

#endif /* ENZO_ENZO_EOS_IDEAL_HPP */
