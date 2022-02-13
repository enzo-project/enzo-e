// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemannImplNew.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 16 2019
/// @brief    [\ref Enzo] Implementation of RiemannImpl, which is a class
/// template that can be specialized to implement various Riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_IMPL_HPP
#define ENZO_ENZO_RIEMANN_IMPL_HPP

#include <pup_stl.h>
#include <type_traits>
#include <functional> // used to check function signature

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

template<typename EOSStructT>
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

  // EOS Struct Object
  const EOSStructT eos;

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
class EnzoRiemannImpl : public EnzoRiemann
{
  /// @class    EnzoRiemannImpl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides implementation of approximate Riemann
  ///           Solvers
  ///
  /// @tparam ImplFunctor The functor used to specialize `EnzoRiemannImpl`. The
  ///     functor must provide a public member type called `LUT`, which is a
  ///     specialization of `EnzoRiemannLUT<InputLUT>`. It must also support
  ///     initialization from KernelConfig and provide a public `operator()`
  ///     method which accepts 3 integer indices.
  ///
  /// EnzoRiemannImpl factors out the repeated code between different
  /// approximate Riemann Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF &
  /// Roe solvers).
  ///
  /// To define a new RiemannSolver using `EnzoRiemann`:
  ///   1. Define a new `ImplFunctor` (e.g. `HLLKernel`).
  ///   2. It might be useful to define an alias name for the specialization of
  ///      `EnzoRiemannImpl` that uses the new `ImplFunctor`
  ///      (e.g. `using EnzoRiemannHLLD = EnzoRiemannImpl<HLLKernel>;`).
  ///   3. Update `EnzoRiemann::construct_riemann` to construct the Riemann
  ///      Solver when the correct name is specified.
  ///   4. Update the documentation with the name of the newly available
  ///      RiemannSolver

  using LUT = typename ImplFunctor::LUT;
  using EOSStructT = typename ImplFunctor::EOSStructT;

  // Check whether ImplFunctor's operator() method has the expected signature
  // and raise an error message if it doesn't:
  static_assert(std::is_convertible<ImplFunctor&&,
                                    std::function<void(int, int, int)>>::value,
                "ImplFunctor must define the method: "
                "void ImplFunctor::operator() (int,int,int)");



public: // interface

  /// Constructor
  ///
  /// @param internal_energy Indicates whether internal_energy is an
  ///     integration quantity.
  EnzoRiemannImpl(const EnzoRiemann::FactoryArgs factory_args,
                  bool internal_energy);

  /// Virtual destructor
  virtual ~EnzoRiemannImpl(){ delete scratch_ptr_; };

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

private: // helper methods

  /// Actually executes the Riemann Solver
  ///
  /// @note
  /// for purposes of getting icc to vectorize code, it seems to be important
  /// that this method's contents are separated from `EnzoRiemannImpl::solve`
  static void solve_(const KernelConfig<EOSStructT> config,
                     const int stale_depth) noexcept;

private: //attributes

  /// expected keys (and key-order) that the `solve` method expects the
  /// `flux_map` argument to have
  std::vector<std::string> integration_quantity_keys_;

  /// expected keys (and key-order) that the `solve` method expects the
  /// `priml_map` and `primr_map` arguments to have (i.e. these are the keys
  /// for the primitives that are required to compute the flux)
  std::vector<std::string> primitive_quantity_keys_;

  /// Tracks whether the internal energy needs to be computed
  bool calculate_internal_energy_flux_;

  /// Holds lazily allocated scratch arrays
  enzo_riemann_utils::ScratchArrays_* scratch_ptr_;
};

//----------------------------------------------------------------------

template <class ImplFunctor>
EnzoRiemannImpl<ImplFunctor>::EnzoRiemannImpl
(const EnzoRiemann::FactoryArgs factory_args,const bool internal_energy)
  : EnzoRiemann(factory_args)
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

  scratch_ptr_ = new enzo_riemann_utils::ScratchArrays_();
}

//----------------------------------------------------------------------

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::solve
(const EnzoEFltArrayMap &prim_map_l, const EnzoEFltArrayMap &prim_map_r,
 EnzoEFltArrayMap &flux_map, const int dim, const EnzoEquationOfState *eos,
 const int stale_depth, const str_vec_t &passive_list,
 const CelloArray<enzo_float,3> * const interface_velocity) const
{

  // Currently just going to assume that we have an Ideal Equation of State
  // (should probably test that...)
  const enzo_float gamma = eos->get_gamma();
  const EOSStructIdeal eos_struct = {gamma};

  static_assert(std::is_same<EOSStructIdeal, EOSStructT>::value,
                "We currently assume that all EOSs are ideal");

  // The strategy is to allocate some scratch space for internal_energy_flux &
  // velocity_i_bar_array, even if we don't care about dual-energy in order to
  // avoid branching.
  const bool calculate_internal_energy_flux = calculate_internal_energy_flux_;
  EFlt3DArray internal_energy_flux, velocity_i_bar_array;
  enzo_riemann_utils::prep_dual_energy_arrays_(calculate_internal_energy_flux,
                                               flux_map, interface_velocity,
                                               scratch_ptr_,
                                               internal_energy_flux,
                                               velocity_i_bar_array);

#ifdef RIEMANN_DEBUG
  check_key_order_(prim_map_l, true, passive_list);
  check_key_order_(prim_map_r, true, passive_list);
  check_key_order_(flux_map, false, passive_list);
#endif

  const KernelConfig<EOSStructT> config = {dim,
                                           flux_map.get_backing_array(),
                                           prim_map_l.get_backing_array(),
                                           prim_map_r.get_backing_array(),
                                           eos_struct,
                                           internal_energy_flux,
                                           velocity_i_bar_array};

  solve_(config, stale_depth);

  enzo_riemann_utils::solve_passive_advection(prim_map_l, prim_map_r, flux_map,
                                              flux_map.at("density"),
                                              stale_depth, passive_list);

  if (LUT::has_bfields){
    // If Dedner Fluxes are required, they might get handled here
  }

}

//----------------------------------------------------------------------

template <class ImplFunctor>
void EnzoRiemannImpl<ImplFunctor>::solve_
(const KernelConfig<EOSStructT> config, const int stale_depth)
  noexcept
{
  const ImplFunctor kernel{config};

  // preload shape of arrays (to inform compiler they won't change)
  const int mz = config.flux_arr.shape(1);
  const int my = config.flux_arr.shape(2);
  const int mx = config.flux_arr.shape(3);

  // compute the flux at all non-stale cell interfaces
  for (int iz = stale_depth; iz < mz - stale_depth; iz++) {
    for (int iy = stale_depth; iy < my - stale_depth; iy++) {
      #pragma omp simd
      for (int ix = stale_depth; ix < mx - stale_depth; ix++) {
        kernel(iz,iy,ix);
      }
    }
  }
}

#endif /* ENZO_ENZO_RIEMANN_IMPL_HPP */
