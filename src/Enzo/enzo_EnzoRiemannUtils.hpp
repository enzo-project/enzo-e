// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannUtils.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 15 2020
/// @brief    [\ref Enzo] Declaration and Implementation of assorted utilities
/// in the enzo_riemann_utils namesapce that are associated with the
/// implementation of Riemann Solvers.
///
/// Many of these functions could conceivably be methods of EnzoRiemannImpl,
/// but doing so would cause the compiler to generate unnecessary duplicate
/// versions of of template functions for the different ImplFunctors that use
/// the same underlying LUT. It might be worth considering making these
/// functions into static methods of EnzoRiemannLUT.

#ifndef ENZO_ENZO_RIEMANN_UTILS_HPP
#define ENZO_ENZO_RIEMANN_UTILS_HPP

//----------------------------------------------------------------------

struct EOSStructIdeal{
  // In the future it might make sense to use this under the hood of the
  // EnzoEOSIdeal class.
  //
  // It might also make sense to make this into POD struct that get's passed to
  // an Impl Functor's operator() method. If the structure is simple enough,
  // the overhead of constructing/destroying can be optimized out

  static inline enzo_float specific_eint(const enzo_float density,
                                         const enzo_float pressure,
                                         const enzo_float gamma) noexcept
  { return pressure / ( (gamma - 1.0) * density); }

  static inline enzo_float eint_dens(const enzo_float density, const enzo_float pressure,
                                     const enzo_float gamma) noexcept
  { return pressure / (gamma - 1.0); }

  static inline enzo_float sound_speed(const enzo_float density, const enzo_float pressure,
                                       const enzo_float gamma) noexcept
  { return std::sqrt(gamma * pressure / density); }

};

//----------------------------------------------------------------------

namespace enzo_riemann_utils{

  template <class LUT>
  inline enzo_float mag_pressure(const lutarray<LUT> prim) noexcept
  {
    enzo_float bi = (LUT::bfield_i >= 0) ? prim[LUT::bfield_i] : 0;
    enzo_float bj = (LUT::bfield_j >= 0) ? prim[LUT::bfield_j] : 0;
    enzo_float bk = (LUT::bfield_k >= 0) ? prim[LUT::bfield_k] : 0;
    return 0.5 * (bi*bi + bj*bj + bk *bk);
  }

  //----------------------------------------------------------------------

  /// Computes the conserved counterpart for every primitive.
  ///
  /// Primitives are generally categorized as conserved or specific. A
  /// "conserved" primitive is directly copied, and a "specific" primitive
  /// is multiplied by the density. The main exception is pressure, which is
  /// stored in prim[LUT::total_energy].
  ///
  /// This is currently assumed to only be implemented for a non-barotropic EOS
  /// with a calorically perfect EOS
  ///
  /// We might want to consolidate this with active_fluxes
  template <class LUT>
  inline lutarray<LUT> compute_conserved(const lutarray<LUT> prim,
                                         const enzo_float gamma) noexcept
  {
    lutarray<LUT> cons;

    // the values held in prim at index 0 up to (but not including)
    // LUT::specific_start are already conserved quantities.
    for (std::size_t i = 0; i < LUT::specific_start; i++) {
      cons[i] = prim[i];
    }

    // the remainder of the values are specific quantities
    for (std::size_t i = LUT::specific_start; i < LUT::num_entries; i++) {
      cons[i] = prim[i] * prim[LUT::density];
    }

    if (LUT::total_energy >= 0) { // overwrite the total energy index
      enzo_float density = prim[LUT::density];
      enzo_float pressure = prim[LUT::total_energy];
      enzo_float internal_edens = EOSStructIdeal::eint_dens(density, pressure,
                                                            gamma);

      const enzo_float vi = prim[LUT::velocity_i];
      const enzo_float vj = prim[LUT::velocity_j];
      const enzo_float vk = prim[LUT::velocity_k];
      const enzo_float kinetic_edens = 0.5 * density * (vi*vi + vj*vj + vk*vk);

      enzo_float magnetic_edens = mag_pressure<LUT>(prim);

      cons[LUT::total_energy] = internal_edens + kinetic_edens + magnetic_edens;
    }
    
    return cons;
  }

  //----------------------------------------------------------------------

  template <class LUT>
  inline enzo_float fast_magnetosonic_speed(const lutarray<LUT> prim_vals,
                                            enzo_float pressure,
                                            enzo_float gamma) noexcept
  {
    enzo_float bi = (LUT::bfield_i >= 0) ? prim_vals[LUT::bfield_i] : 0;
    enzo_float bj = (LUT::bfield_j >= 0) ? prim_vals[LUT::bfield_j] : 0;
    enzo_float bk = (LUT::bfield_k >= 0) ? prim_vals[LUT::bfield_k] : 0;

    // TODO: optimize calc of cs2 to omit sqrt and pow
    //       can also skip the calculation of B2 by checking if
    //       LUT::bfield_i, LUT::bfield_j, LUT::bfield_k are all negative
    const enzo_float cs = EOSStructIdeal::sound_speed(prim_vals[LUT::density],
                                                pressure, gamma);
    const enzo_float cs2 = std::pow(cs,2);
    const enzo_float B2 = (bi*bi + bj*bj + bk *bk);
    if (B2 == 0){
      return std::sqrt(cs2);
    }
    enzo_float va2 = B2/prim_vals[LUT::density];
    // TODO: replace va2 * cos2 with va2_cos2 = bi*bi/prim_vals[LUT::density]
    enzo_float cos2 = bi*bi / B2;
    return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
                                            4.*cs2*va2*cos2)));
  }

  //----------------------------------------------------------------------

  /// This function should be called when we to fix cos2 to some predetermined
  /// value
  template <class LUT>
  inline enzo_float fast_magnetosonic_speed(const lutarray<LUT> prim_vals,
                                            enzo_float pressure,
                                            enzo_float gamma,
                                            enzo_float cos2) noexcept
  {
    enzo_float bi = (LUT::bfield_i >= 0) ? prim_vals[LUT::bfield_i] : 0;
    enzo_float bj = (LUT::bfield_j >= 0) ? prim_vals[LUT::bfield_j] : 0;
    enzo_float bk = (LUT::bfield_k >= 0) ? prim_vals[LUT::bfield_k] : 0;

    // TODO: optimize calc of cs2 to omit sqrt and pow
    //       can also skip the calculation of B2 by checking if
    //       LUT::bfield_i, LUT::bfield_j, LUT::bfield_k are all negative
    const enzo_float cs = EOSStructIdeal::sound_speed(prim_vals[LUT::density],
                                                pressure, gamma);
    const enzo_float cs2 = std::pow(cs,2);
    const enzo_float B2 = (bi*bi + bj*bj + bk *bk);
    if (B2 == 0){
      return std::sqrt(cs2);
    }
    enzo_float va2 = B2/prim_vals[LUT::density];
    return std::sqrt(0.5*(va2+cs2+std::sqrt(std::pow(cs2+va2,2) -
                                            4.*cs2*va2*cos2)));
  }

  //----------------------------------------------------------------------

  /// computes fluxes for the basic mhd conserved quantities - density,
  /// momentum, energy, magnetic fields
  template <class LUT>
  inline lutarray<LUT> active_fluxes(const lutarray<LUT> prim,
                                     const lutarray<LUT> cons,
                                     enzo_float pressure) noexcept
  {
    lutarray<LUT> fluxes;
    enzo_float vi, vj, vk, p, Bi, Bj, Bk, etot, magnetic_pressure;
    vi = prim[LUT::velocity_i];
    vj = prim[LUT::velocity_j];
    vk = prim[LUT::velocity_k];

    Bi = (LUT::bfield_i >= 0) ? prim[LUT::bfield_i] : 0;
    Bj = (LUT::bfield_j >= 0) ? prim[LUT::bfield_j] : 0;
    Bk = (LUT::bfield_k >= 0) ? prim[LUT::bfield_k] : 0;
    etot =  (LUT::total_energy >= 0) ? cons[LUT::total_energy] : 0;
    
    p = pressure;

    magnetic_pressure = mag_pressure<LUT>(prim);

    // Compute Fluxes
    enzo_float mom_i = cons[LUT::velocity_i];
    fluxes[LUT::density] = mom_i;

    // Fluxes for Mx, My, Mz
    fluxes[LUT::velocity_i] = mom_i*vi - Bi*Bi + p + magnetic_pressure;
    fluxes[LUT::velocity_j] = mom_i*vj - Bj*Bi;
    fluxes[LUT::velocity_k] = mom_i*vk - Bk*Bi;

    // Flux for etot
    fluxes[LUT::total_energy] = ((etot + p + magnetic_pressure)*vi
                                 - (Bi*vi + Bj*vj + Bk*vk)*Bi);

    // Fluxes for Bi,Bj,Bk
    if (LUT::bfield_i >= 0) { fluxes[LUT::bfield_i] = 0; }
    if (LUT::bfield_j >= 0) { fluxes[LUT::bfield_j] = Bj*vi - Bi*vj; }
    if (LUT::bfield_k >= 0) { fluxes[LUT::bfield_k] = Bk*vi - Bi*vk; }

    return fluxes;
  }

  //----------------------------------------------------------------------

  inline std::string parse_mem_name_(std::string member_name,
                                     EnzoPermutedCoordinates coord)
  {
    char component;

    // this helper function is ONLY invoked in cases where we know member_name
    // is guaranteed to correspond to an actively advected quantity (so we
    // don't need to check the return value of the following function)
    EnzoCenteredFieldRegistry::get_actively_advected_quantity_name
      (member_name, true, &component);

    // return immediately, if not a vector component
    if (component == '\0') { return member_name; }

    std::string out = member_name.substr(0, member_name.length() - 1);
    char suffixes[3] {'x','y','z'};
    if (component == 'i'){
      out.push_back(suffixes[coord.i_axis()]);
    } else if (component == 'j'){
      out.push_back(suffixes[coord.j_axis()]);
    } else if (component == 'k'){
      out.push_back(suffixes[coord.k_axis()]);
    } else {
      ERROR("enzo_riemann_utils::parse_mem_name_",
	    "branch should be unreachable");
    }
    return out;
  }

  //----------------------------------------------------------------------

  /// Returns the keys corresponding to each entry in the LUT, (in the
  /// order of the lookup table). For components of vector quantities, this
  /// maps the i, j, and k components to the x, y, and z components
  ///
  /// @param prim When true, the function returns the keys for the primitive
  ///   quantities. Otherwise, it returns the keys for the integration
  ///   quantities.
  template<class LUT>
  inline std::vector<std::string> get_quantity_keys(const bool prim) noexcept
  {
    // initialize out as a vector of empty strings
    std::vector<std::string> out(LUT::num_entries);

    EnzoPermutedCoordinates coord(0); // map i,j,k to x,y,z
    // define a lambda function to execute for every member of lut. For each
    // member in lut, its passed: 1. the member's name
    //                            2. the associated index
    auto fn = [coord, prim, &out](const std::string& name, const int index)
      {
        if (index != -1){
          if (prim && (index == LUT::total_energy)){
            out[index] = "pressure";
          } else {
            out[index] = parse_mem_name_(name, coord);
          }
        }
      };

    LUT::for_each_entry(fn);

    return out;
  }

  //----------------------------------------------------------------------

  /// Constructs an array of instances of EFlt3DArray this is something of a
  /// stopgap solution until we fully support implemntation of EFlt3DArray
  /// with 4D CelloArrays
  template<class LUT>
  inline std::array<CelloArray<const enzo_float, 3>, LUT::num_entries>
  array_from_map(const EnzoEFltArrayMap& map) noexcept
  {
    std::array<CelloArray<const enzo_float, 3>, LUT::num_entries> out;
    for (std::size_t i = 0; i < LUT::num_entries; i++) { out[i] = map[i]; }
    return out;
  }

  template<class LUT>
  inline std::array<CelloArray<enzo_float, 3>,LUT::num_entries> array_from_map
  (EnzoEFltArrayMap& map) noexcept
  {
    std::array<CelloArray<enzo_float, 3>, LUT::num_entries> out;
    for (std::size_t i = 0; i < LUT::num_entries; i++) { out[i] = map[i]; }
    return out;
  }

  //----------------------------------------------------------------------

  /// @def      LUT_TRANSFER_T_SCALAR
  /// @brief    Part of the LUT_TRANSFER group of macros that are used to copy
  ///           scalar actively advected values between a
  ///           std::array<enzo_float, LUT::num_entries> and a given
  ///           location in a list of 3D external arrays. These macros simply
  ///           ignore the vector quantities (since that involves some level of
  ///           permutation)
  ///
  /// The macro copies values from `src` to `dest`. When `src` (`dest`) is the
  /// list of external arrays, `SRCSUF` (`DESTSUF`) should be the parentheses
  /// enclosed indices of the arrays while `DESTSUF` (`SRCSUF`) should be passed
  /// an empty argument.
  #define LUT_TRANSFER_T_SCALAR(LUT, src, dest, SRCSUF, DESTSUF, name)        \
    if (LUT::name != -1){ dest[LUT::name] DESTSUF = src[LUT::name] SRCSUF; }
  #define LUT_TRANSFER_T_VECTOR(LUT, src, dest, SRCSUF, DESTSUF, name) /*...*/
  #define LUT_TRANSFER_F_SCALAR(LUT, src, dest, SRCSUF, DESTSUF, name) /*...*/
  #define LUT_TRANSFER_F_VECTOR(LUT, src, dest, SRCSUF, DESTSUF, name) /*...*/

  //----------------------------------------------------------------------

  /// transfers actively advected scalar quantities values from the specified
  /// correct location in `external` to a 1D output array.
  template<class LUT>
  void transfer_scalars_to_lutarray
  (const int dim, const int iz, const int iy, const int ix,
   const std::array<CelloArray<const enzo_float, 3>,LUT::num_entries>& external,
   lutarray<LUT> &dest)
    noexcept
  {
    #define ENTRY(name, math_type, category, if_advection)                    \
      LUT_TRANSFER_##if_advection##_##math_type (EnzoRiemannLUT<LUT>,         \
                                                 external, dest, (iz,iy,ix), ,\
                                                 name);
    FIELD_TABLE
    #undef ENTRY
  }

  //----------------------------------------------------------------------

  /// transfers actively advected scalar quantities values from `src` to the
  /// correct location in `external`.
  template<class LUT>
  inline void transfer_scalars_from_lutarray
  (const int dim, const int iz, const int iy, const int ix,
   const std::array<CelloArray<enzo_float, 3>, LUT::num_entries>& external,
   const lutarray<LUT> src) noexcept
  {
    #define ENTRY(name, math_type, category, if_advection)                    \
      LUT_TRANSFER_##if_advection##_##math_type (EnzoRiemannLUT<LUT>,         \
                                                 src, external, , (iz,iy,ix), \
                                                 name);
    FIELD_TABLE
    #undef ENTRY
  }

}

#endif /* ENZO_ENZO_RIEMANN_UTILS_HPP */
