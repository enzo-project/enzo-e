// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannImplNew.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 16 2019
/// @brief    [\ref Enzo] Implementation of RiemannImpl, which is a class
/// template that can be specialized to implement various Riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_IMPL2_HPP
#define ENZO_ENZO_RIEMANN_IMPL2_HPP

#include <pup_stl.h>
#include <cstdint> // used to check that static methods are defined

// defining RIEMANN_DEBUG adds some extra error checking, useful for debugging
//#define RIEMANN_DEBUG

//----------------------------------------------------------------------

struct KernelConfig{
  /// @class    KernelConfig
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Stores the configuration options that would be
  /// passed to a Riemann Solver Function
  ///
  /// to help facilitate optimizations, avoid introducing methods and make all
  /// member variables ``const``

  /// dimension along which the flux is being computed
  /// 0 -> x, 1 -> y, 2 -> z 
  const int dim;

  /// arrays to store the computed fluxes
  const CelloArray<enzo_float,4> flux_arr;
  /// arrays storing the left and right reconstructed primitives
  const CelloArray<const enzo_float,4> prim_arr_l;
  const CelloArray<const enzo_float,4> prim_arr_r;

  // EOS related parameters:

  /// specifies the adiabatic index (if we support barotropic EOSs, we may need
  /// to support isothermal sound-speed in the future)
  const enzo_float gamma;

  // Dual-Energy Related Parameters:

  /// array to store the computed fluxes for the specific internal energy (this
  /// may alias an array in fluxes)
  const CelloArray<enzo_float,3> internal_energy_flux_arr;
  /// array to hold the computed component of the velocity at the cell
  /// interfaces along `dim`
  const CelloArray<enzo_float,3> velocity_i_bar_arr;

};

//----------------------------------------------------------------------

template <class ImplFunctor>
class EnzoRiemannImplNew : public EnzoRiemann
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

  //static_assert(std::is_default_constructible<ImplFunctor>::value,
  //		"ImplFunctor is not default constructable");

  // Check whether ImplFunctor's operator() method has the expected signature
  // and raise an error message if it doesn't:
  //   - first, define a struct (has_expected_functor_sig_) to check if
  //     ImplFunctor's operator() method has the expected signature. This needs
  //     to happen in the current scope because the signature depends on LUT
  //DEFINE_HAS_INSTANCE_METHOD(has_expected_functor_sig_, operator(),
  //                           riemann_function_call_signature<LUT>);
  //static_assert(has_expected_functor_sig_<ImplFunctor>::value,
  //		"ImplFunctor's operator() method doesn't have the correct "
  //		"function signature");


public: // interface

  /// Constructor
  ///
  /// @param internal_energy Indicates whether internal_energy is an
  ///     integration quantity.
  EnzoRiemannImplNew(const EnzoRiemann::FactoryArgs factory_args,
                   bool internal_energy);

  /// Virtual destructor
  virtual ~EnzoRiemannImplNew(){ };

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
EnzoRiemannImplNew<ImplFunctor>::EnzoRiemannImplNew
(const EnzoRiemann::FactoryArgs factory_args,const bool internal_energy)
  : EnzoRiemann(factory_args)
{
  integration_quantity_keys_ =
    enzo_riemann_utils::get_quantity_keys<LUT>(false);
  primitive_quantity_keys_ = enzo_riemann_utils::get_quantity_keys<LUT>(true);

  for (std::string key : integration_quantity_keys_){
    ASSERT("EnzoRiemannImplNew::EnzoRiemannImplNew",
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
void EnzoRiemannImplNew<ImplFunctor>::solve
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
  ASSERT("EnzoRiemannImplNew::solve", "currently no support for barotropic eos",
	 !barotropic);

  // preload shape of arrays (to inform compiler they won't change)
  const int mz = flux_map.array_shape(0);
  const int my = flux_map.array_shape(1);
  const int mx = flux_map.array_shape(2);

  // Check if the internal energy flux must be computed and then initialize the
  // appropriate array
  EFlt3DArray internal_energy_flux, velocity_i_bar_array;
  if ((calculate_internal_energy_flux_) && (interface_velocity == nullptr)){
    ERROR("EnzoRiemannImplNew::solve",
          "interface_velocity is expected to be non-NULL when computing the "
          "internal energy flux");
  } else if (calculate_internal_energy_flux_) {
    internal_energy_flux = flux_map.at("internal_energy");
    velocity_i_bar_array = *interface_velocity;
  } else {
    // TODO: preallocate scratch space
    internal_energy_flux = EFlt3DArray(mz,my,mx);
    velocity_i_bar_array = EFlt3DArray(mz,my,mx);
  }

#ifdef RIEMANN_DEBUG
  check_key_order_(prim_map_l, true, passive_list);
  check_key_order_(prim_map_r, true, passive_list);
  check_key_order_(flux_map, false, passive_list);
#endif

  const KernelConfig config{dim,
                            flux_map.get_backing_array(),
                            prim_map_l.get_backing_array(),
                            prim_map_r.get_backing_array(),
                            gamma,
                            internal_energy_flux,
                            velocity_i_bar_array};
  const ImplFunctor kernel{config};

  // compute the flux at all non-stale cell interfaces
  const int sd = stale_depth;
  for (int iz = sd; iz < mz - sd; iz++) {
    for (int iy = sd; iy < my - sd; iy++) {
      #pragma omp simd
      for (int ix = sd; ix < mx - sd; ix++) {
        kernel(iz,iy,ix);
      }
    }
  }

  enzo_riemann_utils::solve_passive_advection(prim_map_l, prim_map_r, flux_map,
                                              flux_map.at("density"),
                                              stale_depth, passive_list);

  if (LUT::has_bfields){
    // If Dedner Fluxes are required, they might get handled here
  }

}

#endif /* ENZO_ENZO_RIEMANN_IMPL2_HPP */
