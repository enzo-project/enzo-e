// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannImpl.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 16 2019
/// @brief    [\ref Enzo] Implementation of RiemannImpl, which is a class
/// template that can be specialized to implement various Riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_IMPL_HPP
#define ENZO_ENZO_RIEMANN_IMPL_HPP

#include <pup_stl.h>
#include <cstdint> // used to check that static methods are defined

// defining RIEMANN_DEBUG adds some extra error checking, useful for debugging
//#define RIEMANN_DEBUG

//----------------------------------------------------------------------

struct HydroLUT {
  /// @class    HydroLUT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates a compile-time LUT for pure
  ///           hydrodynamics that is used for implementing Riemann Solvers

  enum vals { density=0,
	      velocity_i,
	      velocity_j,
	      velocity_k,
	      total_energy,
	      num_entries};

  // in the future, automatically calculate the following
  static const std::size_t specific_start = 1;
};

//----------------------------------------------------------------------

struct MHDLUT{
  /// @class    MHDLUT
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates a compile-time LUT for MHD that is
  ///           that is to be used for implementing Riemann Solvers
  enum vals { density=0,
	      bfield_i,
	      bfield_j,
	      bfield_k,
	      velocity_i,
	      velocity_j,
	      velocity_k,
	      total_energy,
	      num_entries};
  // in the future, automatically calculate the following
  static const std::size_t specific_start = 4;
};

//----------------------------------------------------------------------

/// @typedef riemann_function_call_signature
/// @brief   This is the function call signature of the function call operator,
///          `operator()`, that is expected for an `ImplFunctor` used to
///          implement a Riemann Solver by specializing `EnzoRiemannImpl`.
///
/// An `ImplFunctor`'s `operator()` method should be declared as:
///
/// @code
///     lutarray<LUT> operator()
///       (const lutarray<LUT> flux_l, const lutarray<LUT> flux_r,
///        const lutarray<LUT> prim_l, const lutarray<LUT> prim_r,
///        const lutarray<LUT> cons_l, const lutarray<LUT> cons_r,
///        enzo_float pressure_l, enzo_float pressure_r, bool barotropic_eos,
///        enzo_float gamma, enzo_float isothermal_cs, enzo_float &vi_bar)
/// @endcode
///
/// This function should computes the Riemann Flux at a given cell-interface
/// for the set of actively advected quantities designated by the `LUT` type
/// publicly defined within the scope of `ImplFunctor` (that can be acessed with
/// `ImplFunctor::LUT`).
///
/// The arguments, `flux_l`, `flux_r`, `cons_l`, and `cons_r` are constant
/// arrays that store values associated with each actively advected integration
/// quantity, at the indices designated by `ImplFunctor::LUT`. For each
/// quantity, `flux_l`/`flux_r` hold the associated fluxes and `cons_l`/`cons_r`
/// hold the associated values (in conserved form). The arguments `prim_l` and
/// `prim_r` are similar, but hold the associated primitive values instead. The
/// indices are also deignated by `ImplFunctor::LUT`. Note that if
/// `ImplFunctor::LUT::total_energy` points to a valid index, then
/// `prim_l`/`prim_r` store pressure at that location.
///
/// The left and right reconstructed pressure values are passed as `pressure_l`
/// and `pressure_r`. `barotropic_eos` indicates whether the fluid equation of
/// state is barotropic. If `true`, then `isothermal_cs` is expected to be
/// non-zero and if `false`, then `gamma` is expected to be positive. `vi_bar`
/// is the estimate of `velocity_i` at the cell-interface to be used in the
/// internal energy source term (for the dual energy formalism).
///
/// The function is expected to return the computed fluxes for each actively
/// advected quantity in an array with mapping set by `ImplFunctor::LUT`
template<class LUT>
using riemann_function_call_signature =
  lutarray<LUT>(*)(const lutarray<LUT> flux_l, const lutarray<LUT> flux_r,
                   const lutarray<LUT> prim_l, const lutarray<LUT> prim_r,
                   const lutarray<LUT> cons_l, const lutarray<LUT> cons_r,
                   enzo_float pressure_l, enzo_float pressure_r,
                   bool barotropic_eos, enzo_float gamma,
                   enzo_float isothermal_cs, enzo_float &vi_bar);

//----------------------------------------------------------------------

template <class ImplFunctor>
class EnzoRiemannImpl : public EnzoRiemann
{
  /// @class    EnzoRiemannImpl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides implementation of approximate Riemann
  ///           Solvers
  ///
  /// @tparam ImplFunctor The functor used to specialize `EnzoRiemannImpl`. The
  ///     functor must provide a public member type called `LUT`, which is a
  ///     specialization of `EnzoRiemannLUT<InputLUT>`. It must also be default
  ///     constructible and provide a public `operator()` method. For details
  ///     about the method's expected signature, see the documentation for
  ///     `riemann_function_call_signature`.
  ///
  /// EnzoRiemannImpl factors out the repeated code between different
  /// approximate Riemann Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF &
  /// Roe solvers).
  ///
  /// To define a new RiemannSolver using `EnzoRiemann`:
  ///   1. Define a new `ImplFunctor` (e.g. `HLLDImpl`).
  ///   2. It might be useful to define an alias name for the specialization of
  ///      `EnzoRiemannImpl` that uses the new `ImplFunctor`
  ///      (e.g. `using EnzoRiemannHLLD = EnzoRiemannImpl<HLLDImpl>;`).
  ///   3. Then add the particlular specialization to enzo.CI (e.g. add the
  ///      line: `PUPable EnzoRiemannImpl<HLLDImpl>;`)
  ///   4. Update `EnzoRiemann::construct_riemann` to construct the Riemann
  ///      Solver when the correct name is specified.
  ///   5. Update the documentation with the name of the newly available
  ///      RiemannSolver

  using LUT = typename ImplFunctor::LUT;

  static_assert(std::is_default_constructible<ImplFunctor>::value,
		"ImplFunctor is not default constructable");

  // Check whether ImplFunctor's operator() method has the expected signature
  // and raise an error message if it doesn't:
  //   - first, define a struct (has_expected_functor_sig_) to check if
  //     ImplFunctor's operator() method has the expected signature. This needs
  //     to happen in the current scope because the signature depends on LUT
  DEFINE_HAS_INSTANCE_METHOD(has_expected_functor_sig_, operator(),
                             riemann_function_call_signature<LUT>);
  static_assert(has_expected_functor_sig_<ImplFunctor>::value,
		"ImplFunctor's operator() method doesn't have the correct "
		"function signature");


public: // interface

  /// Constructor
  ///
  /// @param internal_energy Indicates whether internal_energy is an
  ///     integration quantity.
  EnzoRiemannImpl(bool internal_energy);

  /// Virtual destructor
  virtual ~EnzoRiemannImpl(){ };

  /// CHARM++ PUP::able declaration
  PUPable_decl_template(EnzoRiemannImpl<ImplFunctor>);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemannImpl (CkMigrateMessage *m)
    : EnzoRiemann(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  void solve (const EnzoEFltArrayMap &prim_map_l,
	      const EnzoEFltArrayMap &prim_map_r,
              EnzoEFltArrayMap &flux_map, const int dim,
	      const EnzoEquationOfState *eos, const int stale_depth,
	      const str_vec_t &passive_list,
              const CelloArray<enzo_float,3> * const interface_velocity) const;

  const std::vector<std::string> integration_quantity_keys() const noexcept
  { return integration_quantity_keys_; }

  const std::vector<std::string> primitive_quantity_keys() const noexcept
  { return primitive_quantity_keys_; }

protected : //methods
  
  /// Computes the fluxes for the passively advected quantites.
  void solve_passive_advection_(const EnzoEFltArrayMap &prim_map_l,
                                const EnzoEFltArrayMap &prim_map_r,
                                EnzoEFltArrayMap &flux_map,
				const EFlt3DArray &density_flux,
                                const int stale_depth,
                                const str_vec_t &passive_list) const throw();

  /// debugging method (checks that the order of keys matches expectations
  void check_key_order_(const EnzoEFltArrayMap &map, bool prim,
			const str_vec_t &passive_list) const noexcept;

protected: //attributes

  /// expected keys (and key-order) that the `solve` method expects the
  /// `flux_map` argument to have
  std::vector<std::string> integration_quantity_keys_;

  /// expected keys (and key-order) that the `solve` method expects the
  /// `priml_map` and `primr_map` arguments to have (i.e. these are the keys
  /// for the primitives that are required to compute the flux)
  std::vector<std::string> primitive_quantity_keys_;

  /// Tracks whether the internal energy needs to be computed
  bool calculate_internal_energy_flux_;
};

//----------------------------------------------------------------------

template <class ImplFunctor>
EnzoRiemannImpl<ImplFunctor>::EnzoRiemannImpl(const bool internal_energy)
  : EnzoRiemann()
{
  integration_quantity_keys_ =
    enzo_riemann_utils::get_quantity_keys<LUT>(false);
  primitive_quantity_keys_ = enzo_riemann_utils::get_quantity_keys<LUT>(true);

  for (std::string key : integration_quantity_keys_){
    ASSERT("EnzoRiemannImpl::EnzoRiemannImpl",
	   "No support for a LUT directly containing \"internal_energy\"",
	   key != "internal_energy");
  }

  calculate_internal_energy_flux_ = internal_energy;
  if (calculate_internal_energy_flux_) {
    integration_quantity_keys_.push_back("internal_energy");
  }
}

//----------------------------------------------------------------------

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::pup (PUP::er &p)
{
  EnzoRiemann::pup(p);

  p|integration_quantity_keys_;
  p|primitive_quantity_keys_;
  p|calculate_internal_energy_flux_;
}

//----------------------------------------------------------------------

/// computes the flux of a passive scalar at a given cell interface.
///
/// @param left,right The passive scalars (as mass fractions) on the left and
///    right sides of the cell interface
/// @param density_flux The density flux at the local interface
static inline enzo_float calc_passive_scalar_flux_(const enzo_float left,
                                                   const enzo_float right,
                                                   const enzo_float density_flux){
  // next line is equivalent to: upwind = (density_flux > 0) ? left : right;
  // but is branchless
  enzo_float upwind = (density_flux > 0) * left + (density_flux <= 0) * right;
  return upwind * density_flux;
}

//----------------------------------------------------------------------

static inline enzo_float passive_eint_flux_(const enzo_float density_l,
                                            const enzo_float pressure_l,
                                            const enzo_float density_r,
                                            const enzo_float pressure_r,
                                            const enzo_float gamma,
                                            const enzo_float density_flux){
  enzo_float eint_l = EOSStructIdeal::specific_eint(density_l, pressure_l,
                                                    gamma);
  enzo_float eint_r = EOSStructIdeal::specific_eint(density_r, pressure_r,
                                                    gamma);
  return calc_passive_scalar_flux_(eint_l, eint_r, density_flux);
}

//----------------------------------------------------------------------

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::solve
(const EnzoEFltArrayMap &prim_map_l, const EnzoEFltArrayMap &prim_map_r,
 EnzoEFltArrayMap &flux_map, const int dim, const EnzoEquationOfState *eos,
 const int stale_depth, const str_vec_t &passive_list,
 const CelloArray<enzo_float,3> * const interface_velocity) const
{

  const bool barotropic = eos->is_barotropic();
  // if not barotropic then the following doesn't have to be reasonable
  const enzo_float isothermal_cs = eos->get_isothermal_sound_speed();
  // If barotropic, then the following doesn't have to be a reasonable value
  // (this will have to be tweaked when we introduce species that modify gamma)
  const enzo_float gamma = eos->get_gamma();

  // When barotropic equations of state are eventually introduced, all eos
  // dependencies should be moved up here
  ASSERT("EnzoRiemannImpl::solve", "currently no support for barotropic eos",
	 !barotropic);

  // TODO: Add special handling for barotropic equations of state
  const CelloArray<const enzo_float, 3> pressure_array_l =
    prim_map_l.at("pressure");
  const CelloArray<const enzo_float, 3> pressure_array_r =
    prim_map_r.at("pressure");

  // Check if the internal energy flux must be computed and then initialize the
  // appropriate array
  const bool calculate_internal_energy_flux = calculate_internal_energy_flux_;
  if ((calculate_internal_energy_flux) && (interface_velocity == nullptr)){
    ERROR("EnzoRiemannImpl::solve",
          "interface_velocity is expected to be non-NULL when computing the "
          "internal energy flux");
  }
  // TODO: check if falling back to default constructor actually introduces a
  // branch in the for-loop. (if so consider using scratch space, instead)
  const EFlt3DArray internal_energy_flux = (calculate_internal_energy_flux) ?
    flux_map.at("internal_energy") : EFlt3DArray();
  const EFlt3DArray velocity_i_bar_array = (calculate_internal_energy_flux) ?
    *interface_velocity : EFlt3DArray();

#ifdef RIEMANN_DEBUG
  check_key_order_(prim_map_l, true, passive_list);
  check_key_order_(prim_map_r, true, passive_list);
  check_key_order_(flux_map, false, passive_list);
#endif

  

  const std::array<CelloArray<const enzo_float, 3>, LUT::num_entries>
    wl_arrays = enzo_riemann_utils::array_from_map<LUT>(prim_map_l);
  const std::array<CelloArray<const enzo_float, 3>, LUT::num_entries> wr_arrays
    = enzo_riemann_utils::array_from_map<LUT>(prim_map_r);
  const std::array<CelloArray<enzo_float, 3>, LUT::num_entries> flux_arrays =
    enzo_riemann_utils::array_from_map<LUT>(flux_map);

  // determine mapping between vector components of the external arrays and the
  // values of LUT (by default, the Riemann Solver maps i->x, j->y, k->z
  const int external_velocity_i = dim + LUT::velocity_i;
  const int external_velocity_j = ((dim+1)%3) + LUT::velocity_i;
  const int external_velocity_k = ((dim+2)%3) + LUT::velocity_i;
  const int external_bfield_i = dim + LUT::bfield_i;
  const int external_bfield_j = ((dim+1)%3) + LUT::bfield_i;
  const int external_bfield_k = ((dim+2)%3) + LUT::bfield_i;

  ImplFunctor func;

  // compute the flux at all non-stale cell interfaces
  const int sd = stale_depth;
  for (int iz = sd; iz < flux_arrays[0].shape(0) - sd; iz++) {
    for (int iy = sd; iy < flux_arrays[0].shape(1) - sd; iy++) {
      for (int ix = sd; ix < flux_arrays[0].shape(2) - sd; ix++) {

        // get the local values of the fluid fields
        lutarray<LUT> wl, wr;
        // first copy the scalar-type quantities (e.g. density, total_energy)
        enzo_riemann_utils::transfer_scalars_to_lutarray<LUT>(dim, iz, iy, ix,
                                                              wl_arrays, wl);
        enzo_riemann_utils::transfer_scalars_to_lutarray<LUT>(dim, iz, iy, ix,
                                                              wr_arrays, wr);
        // now, manually copy the vector-type quantities
        wl[LUT::velocity_i] = wl_arrays[external_velocity_i](iz,iy,ix);
        wl[LUT::velocity_j] = wl_arrays[external_velocity_j](iz,iy,ix);
        wl[LUT::velocity_k] = wl_arrays[external_velocity_k](iz,iy,ix);

        wr[LUT::velocity_i] = wr_arrays[external_velocity_i](iz,iy,ix);
        wr[LUT::velocity_j] = wr_arrays[external_velocity_j](iz,iy,ix);
        wr[LUT::velocity_k] = wr_arrays[external_velocity_k](iz,iy,ix);

        if (LUT::has_bfields()){
          wl[LUT::bfield_i] = wl_arrays[external_bfield_i](iz,iy,ix);
          wl[LUT::bfield_j] = wl_arrays[external_bfield_j](iz,iy,ix);
          wl[LUT::bfield_k] = wl_arrays[external_bfield_k](iz,iy,ix);

          wr[LUT::bfield_i] = wr_arrays[external_bfield_i](iz,iy,ix);
          wr[LUT::bfield_j] = wr_arrays[external_bfield_j](iz,iy,ix);
          wr[LUT::bfield_k] = wr_arrays[external_bfield_k](iz,iy,ix);
        }

        // get the left/right pressure
        enzo_float pressure_l = pressure_array_l(iz,iy,ix);
        enzo_float pressure_r = pressure_array_r(iz,iy,ix);

        // get the conserved quantities
        lutarray<LUT> Ul = enzo_riemann_utils::compute_conserved<LUT>(wl,
                                                                      gamma);
        lutarray<LUT> Ur = enzo_riemann_utils::compute_conserved<LUT>(wr,
                                                                      gamma);

        // compute the interface fluxes
        lutarray<LUT> Fl = enzo_riemann_utils::active_fluxes<LUT>(wl, Ul,
                                                                  pressure_l);
        lutarray<LUT> Fr = enzo_riemann_utils::active_fluxes<LUT>(wr, Ur,
                                                                  pressure_r);

        enzo_float interface_velocity_i;
        // Now compute the Riemann Fluxes
        lutarray<LUT> fluxes = func(Fl, Fr, wl, wr, Ul, Ur, pressure_l,
                                    pressure_r, barotropic, gamma,
                                    isothermal_cs, interface_velocity_i);

        // record the Riemann Fluxes
        // first handle-type scalar quantities (e.g. density, total_energy)
        enzo_riemann_utils::transfer_scalars_from_lutarray<LUT>(dim, iz, iy, ix,
                                                                flux_arrays,
                                                                fluxes);
        // now, manually record fluxes for vector-type quantities
        flux_arrays[external_velocity_i](iz,iy,ix) = fluxes[LUT::velocity_i];
        flux_arrays[external_velocity_j](iz,iy,ix) = fluxes[LUT::velocity_j];
        flux_arrays[external_velocity_k](iz,iy,ix) = fluxes[LUT::velocity_k];
        if (LUT::has_bfields()){
          flux_arrays[external_bfield_i](iz,iy,ix) = fluxes[LUT::bfield_i];
          flux_arrays[external_bfield_j](iz,iy,ix) = fluxes[LUT::bfield_j];
          flux_arrays[external_bfield_k](iz,iy,ix) = fluxes[LUT::bfield_k];
        }

        if (calculate_internal_energy_flux){
          velocity_i_bar_array(iz,iy,ix) = interface_velocity_i;
          internal_energy_flux(iz,iy,ix) = passive_eint_flux_
            (wl[LUT::density], pressure_l, wr[LUT::density], pressure_r,
             gamma, fluxes[LUT::density]);
        }
      }
    }
  }

  solve_passive_advection_(prim_map_l, prim_map_r, flux_map,
  			   flux_arrays[LUT::density], stale_depth,
                           passive_list);

  if (LUT::has_bfields()){
    // If Dedner Fluxes are required, they might get handled here
  }
}

//----------------------------------------------------------------------

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::solve_passive_advection_
(const EnzoEFltArrayMap &prim_map_l, const EnzoEFltArrayMap &prim_map_r,
 EnzoEFltArrayMap &flux_map, const EFlt3DArray &density_flux,
 const int stale_depth, const str_vec_t &passive_list) const throw()
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
  const int sd = stale_depth;
  for (int iz = sd; iz < density_flux.shape(0) - sd; iz++) {
    for (int iy = sd; iy < density_flux.shape(1) - sd; iy++) {

      for (int key_ind=0; key_ind<num_keys; key_ind++){
        for (int ix = sd; ix < density_flux.shape(2) - sd; ix++) {
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

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::check_key_order_
(const EnzoEFltArrayMap &map, bool prim, const str_vec_t &passive_list)
  const noexcept
{
  auto concat = [](str_vec_t vec1, const str_vec_t& vec2){
    vec1.insert( vec1.end(), vec2.begin(), vec2.end() );
    return vec1;
  };

  str_vec_t standard_keys =
    prim ? primitive_quantity_keys() : integration_quantity_keys();
  // confirm that the first `standard_keys.size()` keys of map have the same
  // order as `standard_keys`
  map.validate_key_order(standard_keys, true, true);

  // confirm that all keys in passive_list also appear in map
  // (for now, we'll be permissive about the order of these keys)
  for (const auto& key : passive_list){
    if (!map.contains(key)){
      const std::string &name = map.name();
      ERROR2("EnzoRiemannImpl<ImplFunctor>::check_key_order_",
             "The \"%s\" map is missing a \"%s\" passive scalar key",
             name.c_str(), key.c_str());
    }
  }
}

#endif /* ENZO_ENZO_RIEMANN_IMPL_HPP */
