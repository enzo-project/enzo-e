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

namespace{

  /// computes the squared magnitude of a 3D vector
  ///
  /// @note
  /// This function uses parentheses to dictate the order of operations.
  /// This is primarily done to take advantage of the `-fprotect-parens` flag
  /// provided by several mainstream c++ compilers (e.g. gcc, clang, icpc).
  /// This flag will honor the order of operations specified by parentheses
  /// even when value-unsafe optimizations for floating-point operations are
  /// enabled (e.g. -ffast-math).
  FORCE_INLINE enzo_float squared_mag_vec3D(enzo_float i,
                                            enzo_float j,
                                            enzo_float k) noexcept
  { return ((i*i) + ((j*j) + (k*k))); }

}

//----------------------------------------------------------------------

struct EOSStructIdeal{
  /// @class    EOSStructIdeal
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the equation of state for an ideal gas
  ///           (a.k.a. a calorically perfect gas).
  ///
  /// The plan is to replace EnzoEOSIdeal with new, lighter-weight machinery
  /// that makes use of structs like this one.
  ///
  /// Instances of this class are expected to typically be declared `const`
  ///
  /// @note
  /// Currently all members are public so that this can act as an aggregate in
  /// order to make constructors/destructors trivial (and cheap!).
  ///
  /// @note
  /// For performance purposes, we may want to consider the inverse of
  /// `(gamma - 1)` as a separate member.

public:
  /// stores the adiabtic index
  enzo_float gamma;

private:

  /// computes the sound speed squared
  FORCE_INLINE enzo_float sound_speed_sq_(const enzo_float density,
                                          const enzo_float pressure) const
    noexcept
  { return gamma * pressure / density; }

public:

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
    const enzo_float B2 = squared_mag_vec3D(bfield_i, bfield_j, bfield_k);
    const enzo_float cs2 = sound_speed_sq_(density, pressure);

    // the following branch is evaluated at compile-time
    if ((fixed_cos2 == 0) | (fixed_cos2 == 1)){
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

//----------------------------------------------------------------------

namespace enzo_riemann_utils{

  /// Computes the magnetic pressure
  ///
  /// This returns 0 when the LUT doesn't include magnetic fields
  template <class LUT>
  inline enzo_float mag_pressure(const lutarray<LUT> prim) noexcept
  {
    enzo_float bi = (LUT::bfield_i >= 0) ? prim[LUT::bfield_i] : 0;
    enzo_float bj = (LUT::bfield_j >= 0) ? prim[LUT::bfield_j] : 0;
    enzo_float bk = (LUT::bfield_k >= 0) ? prim[LUT::bfield_k] : 0;
    return 0.5 * squared_mag_vec3D(bi, bj, bk);
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
                                         const EOSStructIdeal& eos) noexcept
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
      enzo_float internal_edens = eos.eint_dens(density, pressure);

      const enzo_float vi = prim[LUT::velocity_i];
      const enzo_float vj = prim[LUT::velocity_j];
      const enzo_float vk = prim[LUT::velocity_k];
      const enzo_float kinetic_edens
        = 0.5 * density * squared_mag_vec3D(vi, vj, vk);

      enzo_float magnetic_edens = mag_pressure<LUT>(prim);

      cons[LUT::total_energy] = internal_edens + kinetic_edens + magnetic_edens;
    }
    
    return cons;
  }

  //----------------------------------------------------------------------

  /// computes the fast magnetosonic speed
  ///
  /// @tparam LUT the lookup table to use with prim_vals
  /// @tparam fixed_cos2 When set to -1, this has no effect. When set to 0 or
  ///     to 1, this fixes cos2 to that value
  template <class LUT, int fixed_cos2 = -1>
  inline enzo_float fast_magnetosonic_speed(const lutarray<LUT> prim_vals,
                                            enzo_float pressure,
                                            const EOSStructIdeal& eos) noexcept
  {
    if (!LUT::has_bfields){
      return eos.sound_speed(prim_vals[LUT::density], pressure);
    } else {
      enzo_float bi = (LUT::bfield_i >= 0) ? prim_vals[LUT::bfield_i] : 0;
      enzo_float bj = (LUT::bfield_j >= 0) ? prim_vals[LUT::bfield_j] : 0;
      enzo_float bk = (LUT::bfield_k >= 0) ? prim_vals[LUT::bfield_k] : 0;
      return eos.fast_magnetosonic_speed<fixed_cos2>(prim_vals[LUT::density],
                                                     pressure, bi, bj, bk);
    }
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
    enzo_float vi = prim[LUT::velocity_i];
    enzo_float vj = prim[LUT::velocity_j];
    enzo_float vk = prim[LUT::velocity_k];

    enzo_float Bi = (LUT::bfield_i >= 0) ? prim[LUT::bfield_i] : 0;
    enzo_float Bj = (LUT::bfield_j >= 0) ? prim[LUT::bfield_j] : 0;
    enzo_float Bk = (LUT::bfield_k >= 0) ? prim[LUT::bfield_k] : 0;
    enzo_float etot = (LUT::total_energy >= 0) ? cons[LUT::total_energy] : 0;

    enzo_float ptot = pressure + mag_pressure<LUT>(prim);

    // Compute Fluxes
    enzo_float mom_i = cons[LUT::velocity_i];
    fluxes[LUT::density] = mom_i;

    // Fluxes for Mx, My, Mz
    fluxes[LUT::velocity_i] = mom_i*vi - Bi*Bi + ptot;
    fluxes[LUT::velocity_j] = mom_i*vj - Bj*Bi;
    fluxes[LUT::velocity_k] = mom_i*vk - Bk*Bi;

    // Flux for etot
    if (LUT::total_energy >= 0) { 
      fluxes[LUT::total_energy] = ((etot + ptot)*vi
                                   - (Bi*vi + (Bj*vj + Bk*vk))*Bi);
    }

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

  /// computes the flux of a passive scalar at a given cell interface.
  ///
  /// @param left,right The passive scalars (as mass fractions) on the left and
  ///    right sides of the cell interface
  /// @param density_flux The density flux at the local interface
  static inline enzo_float calc_passive_scalar_flux_
  (const enzo_float left, const enzo_float right,
   const enzo_float density_flux) noexcept
  {
    // next line is equivalent to: upwind = (density_flux > 0) ? left : right;
    // but is branchless
    enzo_float upwind = (density_flux > 0) * left + (density_flux <= 0) * right;
    return upwind * density_flux;
  }

  //----------------------------------------------------------------------

  /// computes the flux of the specific internal energy at a given cell
  /// interface while treating the specific internal energy as a passive scalar
  static inline enzo_float passive_eint_flux(const enzo_float density_l,
                                             const enzo_float pressure_l,
                                             const enzo_float density_r,
                                             const enzo_float pressure_r,
                                             const EOSStructIdeal& eos,
                                             const enzo_float density_flux)
    noexcept
  {
    enzo_float eint_l = eos.specific_eint(density_l, pressure_l);
    enzo_float eint_r = eos.specific_eint(density_r, pressure_r);
    return calc_passive_scalar_flux_(eint_l, eint_r, density_flux);
  }

  //----------------------------------------------------------------------

  /// compute the flux of the passively advected scalar quantities
  ///
  /// @param[in]  priml_map,primr_map Maps of arrays holding the left/right
  ///     reconstructed face-centered primitives.
  /// @param[out] flux_map Holds arrays where the calculated fluxes for the
  ///     integration quantities will be stored.
  /// @param[in]  density_flux Array of density fluxes at each location. This
  ///     This probably aliases the array returned by `flux_map["density"]`
  /// @param[in]  stale_depth indicates the current stale_depth.
  /// @param[in]  passive_list A list of keys for passive scalars.
  ///
  /// @note
  /// The arrays in `priml_map`, `primr_map`, and `flux_map` should all have
  /// the same shape as `density_flux`.
  static void solve_passive_advection
  (const EnzoEFltArrayMap &prim_map_l, const EnzoEFltArrayMap &prim_map_r,
   EnzoEFltArrayMap &flux_map,
   const CelloArray<const enzo_float,3> &density_flux,
   const int stale_depth, const str_vec_t &passive_list) noexcept
  {
    const std::size_t num_keys = passive_list.size();
    if (num_keys == 0) {return;}

    // This was essentially transcribed from hydro_rk in Enzo:

    // load array of fields
    CelloArray<const enzo_float, 3> *wl_arrays =
      new CelloArray<const enzo_float, 3>[num_keys];
    CelloArray<const enzo_float, 3> *wr_arrays =
      new CelloArray<const enzo_float, 3>[num_keys];
    EFlt3DArray *flux_arrays = new EFlt3DArray[num_keys];

    for (std::size_t ind=0; ind<num_keys; ind++){
      wl_arrays[ind] = prim_map_l.at(passive_list[ind]);
      wr_arrays[ind] = prim_map_r.at(passive_list[ind]);
      flux_arrays[ind] = flux_map.at(passive_list[ind]);
    }

    // preload shape of arrays (to inform compiler they won't change)
    const int mz = density_flux.shape(0);
    const int my = density_flux.shape(1);
    const int mx = density_flux.shape(2);

    for (int iz = stale_depth; iz < mz - stale_depth; iz++) {
      for (int iy = stale_depth; iy < my - stale_depth; iy++) {

        for (std::size_t key_ind = 0; key_ind < num_keys; key_ind++){
          #pragma omp simd
          for (int ix = stale_depth; ix < mx - stale_depth; ix++) {
            const enzo_float dens_flux = density_flux(iz,iy,ix);
            const enzo_float wl = wl_arrays[key_ind](iz,iy,ix);
            const enzo_float wr = wr_arrays[key_ind](iz,iy,ix);
            flux_arrays[key_ind](iz,iy,ix) =
              calc_passive_scalar_flux_(wl, wr, dens_flux);
          }
        }

      }
    }
    delete[] wl_arrays; delete[] wr_arrays; delete[] flux_arrays;

  }


  //----------------------------------------------------------------------

  class ScratchArrays_{

    /// @class    ScratchArrays_
    /// @ingroup  Enzo
    /// @brief    [\ref Enzo] Manages the lifetime of lazily initialized
    ///           scratch arrays
    ///
    /// These arrays are only used to avoid some branching. They will only ever
    /// hold garbage data

  public:
    ScratchArrays_() = default;

    /// Load scratch-space arrays
    ///
    /// If the scratch-space arrays have not been pre-allocated, they will be
    /// allocated by this method
    void get_arrays(int mx, int my, int mz,
                    CelloArray<enzo_float,3>& internal_energy_flux,
                    CelloArray<enzo_float,3>& velocity_i_bar_array) noexcept
    {
      if (internal_energy_flux_.is_null()){
        // make the scratch space bigger than necessary (since the shape will
        // change slightly every time that they are used)
        //
        // given `dim` we could be more precise about the max shape
        internal_energy_flux_ = EFlt3DArray(mz+1,my+1,mx+1);
        velocity_i_bar_array_ = EFlt3DArray(mz+1,my+1,mx+1);
      }

      ASSERT("ScratchArrays_::get_arrays",
             "preallocated scratch-arrays are too small",
             (mz <= internal_energy_flux_.shape(0)) &
             (my <= internal_energy_flux_.shape(1)) &
             (mx <= internal_energy_flux_.shape(2)));

      internal_energy_flux = internal_energy_flux_.subarray
        (CSlice(0, mz),CSlice(0,my),CSlice(0,mx));
      velocity_i_bar_array = velocity_i_bar_array_.subarray
        (CSlice(0, mz),CSlice(0,my),CSlice(0,mx));
    }

  private: // attributes

    CelloArray<enzo_float,3> internal_energy_flux_;
    CelloArray<enzo_float,3> velocity_i_bar_array_;
  };

  //----------------------------------------------------------------------

  /// dumb helper function for setting up arrays for holding the
  /// internal_energy_flux and velocity_i_bar_array
  ///
  /// To simplify things, we want to always provide arrays for storing this
  /// data, even if we don't technically need it, in order to avoid branching
  static void prep_dual_energy_arrays_
  (bool calculate_internal_energy_flux, EnzoEFltArrayMap &flux_map,
   const CelloArray<enzo_float,3> * const interface_velocity,
   ScratchArrays_ * ptr,
   CelloArray<enzo_float,3>& internal_energy_flux,
   CelloArray<enzo_float,3>& velocity_i_bar_array)
  {
    if ((calculate_internal_energy_flux) && (interface_velocity == nullptr)){
      ERROR("EnzoRiemannImpl2::solve",
            "interface_velocity is expected to be non-NULL when computing the "
            "internal energy flux");
    } else if (calculate_internal_energy_flux) {
      internal_energy_flux = flux_map.at("internal_energy");
      velocity_i_bar_array = *interface_velocity;
    } else {
      int mz = flux_map.array_shape(0);
      int my = flux_map.array_shape(1);
      int mx = flux_map.array_shape(2);

      ptr->get_arrays(mx, my, mz, internal_energy_flux, velocity_i_bar_array);
    }
  }

}

#endif /* ENZO_ENZO_RIEMANN_UTILS_HPP */
