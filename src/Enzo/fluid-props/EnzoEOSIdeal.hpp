// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSIdeal.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     2023-04-01
/// @brief    [\ref Enzo] Declaration of the EnzoEOSIdeal struct

#ifndef ENZO_ENZO_EOS_IDEAL_HPP
#define ENZO_ENZO_EOS_IDEAL_HPP

struct EnzoEOSIdeal {
  /// @class    EnzoEOSIdeal
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the equation of state for an ideal gas
  ///           (a.k.a. a calorically perfect gas).
  ///
  /// Instances of this class are expected to typically be declared `const`
  ///
  /// @note
  /// Currently all members are public so that this can act as an aggregate in
  /// order to make constructors/destructors trivial (and cheap!).
  ///
  /// @note
  /// Many of the members of this function that have been annotated as with
  /// ``FORCE_INLINE``, were only annotated that way after verifying that the
  /// annotation is necessary for facillitating certain optimizations within
  /// the Riemann Solvers
  ///
  /// @note
  /// For performance purposes, we may want to consider the inverse of
  /// `(gamma - 1)` as a separate member. If we take this approach, we may want
  /// to use a factory method to build the object and enforce consistency among
  /// the different fields, while maintaining the struct's current status as an
  /// "aggregate". The signature of such a factory method might resemble:
  ///  ``EnzoEOSIdeal EnzoEOSIdeal::build(enzo_float gamma);``

public: // attributes
  // make sure to keep the associated pup function synchronized

  /// stores the adiabtic index
  enzo_float gamma;

public: // public interface common to all EOS types
  constexpr static const char* name() noexcept { return "ideal"; }

  constexpr static bool is_barotropic() noexcept { return false; }

  /// create a string for debugging purposes that represents the EOS's value
  std::string debug_string() const noexcept {
    char buffer[50] = "";
    sprintf(buffer, "EnzoEOSIdeal{ gamma : %#.16g }", gamma);
    return buffer;
  }

private:

  /// computes the sound speed squared
  FORCE_INLINE enzo_float sound_speed_sq_(const enzo_float density,
                                          const enzo_float pressure) const
    noexcept
  { return gamma * pressure / density; }

public:

  /// static method used to build the EnzoEOSIdeal object
  ///
  /// this is a static method (not a constructor) in order to ensure that this
  /// struct is considered an aggregate
  static EnzoEOSIdeal construct(double gamma) noexcept
  {
    ASSERT("EnzoEOSIdeal::construct", "gamma should exceed 1.0", gamma > 1.0);
    return {gamma};
  }

  /// returns the adiabatic index
  FORCE_INLINE enzo_float get_gamma() const noexcept { return gamma; }

  /// computes the specific internal energy
  FORCE_INLINE enzo_float specific_eint(const enzo_float density,
                                        const enzo_float pressure) const
    noexcept
  { return pressure / ( (gamma - 1.0) * density); }

  /// computes the internal energy density
  FORCE_INLINE enzo_float eint_dens(const enzo_float density,
                                    const enzo_float pressure) const noexcept
  { return pressure / (gamma - 1.0); }

  /// computes the adiabatic sound speed
  FORCE_INLINE enzo_float sound_speed(const enzo_float density,
                                      const enzo_float pressure) const noexcept
  { return std::sqrt(sound_speed_sq_(density, pressure)); }

  /// computes the fast magnetosonic speed
  ///
  /// @tparam fixed_cos2 When set to -1, this has no effect. When set to 0 or
  ///     to 1, this fixes cos2 to that value
  ///
  /// The formula for the fast magnetosonic speed typically (cos(theta))^2
  /// where theta is the angle between the velocity component. However, this
  /// method is geared towards 2 main contexts:
  ///    1) determining the signal speed as part of computing the maximum
  ///       timestep (in that case we force (cos(theta))^2 to be 0).
  ///    2) usage within the Riemann solver. In this context, we either want to
  ///       assume that (cos(theta))^2 is 1 OR we assume that the velocity is
  ///       entirely aligned with the ith dimension.
  ///
  /// @note
  /// This method has been implemented so that it will return the correct
  /// answer when all components of the magnetic field are set to 0.
  template<int fixed_cos2 = -1>
  inline enzo_float fast_magnetosonic_speed(const enzo_float density,
                                            const enzo_float pressure,
                                            const enzo_float bfield_i,
                                            const enzo_float bfield_j,
                                            const enzo_float bfield_k)
    const noexcept
  {
    const enzo_float B2 = enzo_utils::squared_mag_vec3D(bfield_i, bfield_j,
                                                        bfield_k);
    const enzo_float cs2 = sound_speed_sq_(density, pressure);

    // the following branch is evaluated at compile-time
    if (fixed_cos2 == 0) {
      enzo_float va2 = B2/density;
      return std::sqrt(va2+cs2);
    } else if (fixed_cos2 == 1) {
      // TODO: we can simplify the returned value, but that will affect the
      //   exact outputs in some of our tests
      enzo_float va2 = B2/density;
      return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
                                              4.*cs2*va2*fixed_cos2)));
    } else {
      const enzo_float inv_density = 1.0/density;
      const enzo_float va2 = B2 * inv_density;
      const enzo_float va2_cos2 = (bfield_i*bfield_i) * inv_density;
      return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
                                              4.*cs2*va2_cos2)));
    }
  }
};

/// function responsible for PUPing
///
/// We have explicitly chosen not to making this a method of the EOS class
/// to try to avoid issues with copying instances to a GPU
inline void pup(PUP::er &p, EnzoEOSIdeal& eos) noexcept {
  // make sure to keep this synchronized with the attributes EnzoEOSIdeal
  p | eos.gamma;
}

#endif /* ENZO_ENZO_EOS_IDEAL_HPP */
