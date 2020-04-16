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
/// versions of of template functions for the different ImplStructs that use
/// the same underlying LUT. It might be worth considering making these
/// functions into static methods of EnzoRiemannLUTWrapper.

#ifndef ENZO_ENZO_RIEMANN_UTILS_HPP
#define ENZO_ENZO_RIEMANN_UTILS_HPP

//----------------------------------------------------------------------

/// @typedef earray
/// @brief   Specialization of std::array to be used to hold enzo_floats for
///          use with Riemann solvers
template<class lut>
using earray = std::array<enzo_float, lut::NEQ>;

//----------------------------------------------------------------------

namespace enzo_riemann_utils{

  /// Computes the conserved counterpart for every integrable primitive.
  /// Integrable primitives are categorized as conserved, specific, and other
  ///
  /// quantities classified as conserved are copied, quantities classified as
  /// primitive are multiplied by density, and for simplicity, quantities
  /// classified as other are copied (There are no obvious cases where there
  /// should ever be a quanitity classified as other
  ///
  /// This has been factored out EnzoRiemannImpl to reduce the number of
  /// template functions created when the same LUT is reused
  
  template <class LUT>
  earray<LUT> compute_conserved(const earray<LUT> prim) noexcept
  {
    earray<LUT> cons;

    // the values held in prim at index 0 up to (but not including)
    // LUT::specific_start are already conserved quantities.
    for (std::size_t i = 0; i < LUT::specific_start; i++) {
      cons[i] = prim[i];
    }

    // the remainder of the values are specific quantities
    for (std::size_t i = LUT::specific_start; i < LUT::NEQ; i++) {
      cons[i] = prim[i] * prim[LUT::density];
    }
    return cons;
  }

  //----------------------------------------------------------------------

  template <class LUT>
  enzo_float mag_pressure(const earray<LUT> prim) noexcept
  {
    enzo_float bi = (LUT::bfield_i >= 0) ? prim[LUT::bfield_i] : 0;
    enzo_float bj = (LUT::bfield_j >= 0) ? prim[LUT::bfield_j] : 0;
    enzo_float bk = (LUT::bfield_k >= 0) ? prim[LUT::bfield_k] : 0;
    return 0.5 * (bi*bi + bj*bj + bk *bk);
  }

  //----------------------------------------------------------------------

  template <class LUT>
  enzo_float sound_speed(const earray<LUT> prim_vals,
                         enzo_float pressure, enzo_float gamma) noexcept
  { return std::sqrt(gamma * pressure / prim_vals[LUT::density]); }

  //----------------------------------------------------------------------

  template <class LUT>
  enzo_float fast_magnetosonic_speed(const earray<LUT> prim_vals,
                                     enzo_float pressure,
                                     enzo_float gamma) noexcept
  {
    enzo_float bi = (LUT::bfield_i >= 0) ? prim_vals[LUT::bfield_i] : 0;
    enzo_float bj = (LUT::bfield_j >= 0) ? prim_vals[LUT::bfield_j] : 0;
    enzo_float bk = (LUT::bfield_k >= 0) ? prim_vals[LUT::bfield_k] : 0;

    // TODO: optimize calc of cs2 to omit sqrt and pow
    //       can also skip the calculation of B2 by checking if
    //       LUT::bfield_i, LUT::bfield_j, LUT::bfield_k are all negative
    enzo_float cs2 = std::pow(sound_speed<LUT>(prim_vals, pressure, gamma),2);
    enzo_float B2 = (bi*bi + bj*bj + bk *bk);
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
  enzo_float fast_magnetosonic_speed(const earray<LUT> prim_vals,
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
    enzo_float cs2 = std::pow(sound_speed<LUT>(prim_vals, pressure, gamma),2);
    enzo_float B2 = (bi*bi + bj*bj + bk *bk);
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
  earray<LUT> active_fluxes(const earray<LUT> prim, const earray<LUT> cons,
                            enzo_float pressure) noexcept
  {
    earray<LUT> fluxes;
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
    enzo_float pi = cons[LUT::velocity_i];
    fluxes[LUT::density] = pi;

    // Fluxes for Mx, My, Mz
    fluxes[LUT::velocity_i] = pi*vi - Bi*Bi + p + magnetic_pressure;
    fluxes[LUT::velocity_j] = pi*vj - Bj*Bi;
    fluxes[LUT::velocity_k] = pi*vk - Bk*Bi;

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

  inline void parse_mem_name_(std::string member_name,
                              EnzoPermutedCoordinates coord,
                              std::string &group_name, int &group_index)
  {
    group_name = member_name;
    group_index = 0;
    std::size_t length = member_name.length();
    if (length >= 2){
      std::string suffix = member_name.substr(length-2,2);
      if (suffix == std::string("_i")){
        group_index = coord.i_axis();
      } else if (suffix == std::string("_j")){
        group_index = coord.j_axis();
      } else if (suffix == std::string("_k")){
        group_index = coord.k_axis();
      } else {
        return;
      }

      group_name = member_name.substr(0,length-2);
    }
  }

  //----------------------------------------------------------------------

  template<class LUT>
  EFlt3DArray* load_array_of_fields(Block *block, Grouping &grouping,
                                    int dim, int stale_depth) noexcept
  {
    std::size_t nfields = LUT::NEQ;
    EFlt3DArray* arr = new EFlt3DArray[nfields];
    // in the case where we don't have reconstructed values (dim = -1) we assume
    // that the that i-axis is aligned with the x-axis
    EnzoPermutedCoordinates coord( (dim == -1) ? 0 : dim);
    EnzoFieldArrayFactory array_factory(block, stale_depth);

    // define a lambda function to execute for every member of lut
    auto fn = [arr, coord, dim, &array_factory, &grouping](std::string name,
                                                           int index)
      {
        // name is the name of a given member of the lut
        // index is the value associated with the member
        if (index != -1){
          int group_index;
          std::string group_name;
          parse_mem_name_(name, coord, group_name, group_index);

          // Sanity Check:
          std::string quantity_type = (group_name == name) ? "SCALAR":"VECTOR";
          int group_size = grouping.size(group_name);

          int expected_size = (quantity_type == "VECTOR") ? 3 : 1;
          ASSERT3("load_array_of_fields_",
                  "The \"%s\" group holds %d field(s). It should hold %d.",
                  group_name.c_str(), group_size, expected_size,
                  group_size == expected_size);

          if (dim != -1){
            arr[index] = array_factory.reconstructed_field(grouping,
                                                           group_name,
                                                           group_index, dim);
          } else {
            arr[index] = array_factory.from_grouping(grouping, group_name,
                                                     group_index);
          }
        }
      };

    unary_lut_for_each_<LUT>(fn);

    return arr;
  }

};

#endif /* ENZO_ENZO_RIEMANN_UTILS_HPP */
