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

// SPHINX-SNIPPET-KERNELCONFIG-START-INCLUDE
template<typename EOSStructT>
struct KernelConfig{
  /// @class    KernelConfig
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Stores the configuration options, input arrays, and
  /// output arrays used by a Riemann Solver Kernel.
  ///
  /// To help facilitate optimizations, we avoid introducing methods and
  /// declare all member variables to be ``const``
  ///
  /// Indices passed to the trailing 3 axes of ANY array (3D & 4D) held by this
  /// struct are used to specify spatial position of values.
  /// - In practice, it's sufficient to understand that a given set of spatial
  ///   indices, `(iz,iy,ix)`, represent the same spatial location in each array
  /// - For completeness, we note that integer indices map to the center of the
  ///   x-faces, y-faces, & z-faces for `dim` values of 0, 1, & 2, respectively.
  ///   The extents of these trailing spatial axes `(LENz, LENy, LENx)` are:
  ///   - when `dim=0`: `LENz = mz`,     `LENy = my`,     `LENx = mx-1`
  ///   - when `dim=1`: `LENz = mz`,     `LENy = my - 1`, `LENx = mx`
  ///   - when `dim=2`: `LENz = mz - 1`, `LENy = my`,     `LENx = mx`
  ///   where `(mz,my,mx)` gives the number of cells along each axis of a
  ///   cell-centered field (including the ghost zone)

  /// @name ConfigOptions
  /**@{*/
  /// dimension along which the flux is being computed (0 -> x, 1 -> y, 2 -> z)
  const int dim;
  /// EOS Struct Object
  const EOSStructT eos;
  /**@}*/

  /// @name PrimaryArrays
  /// `prim_arr_l` and `prim_arr_r` provide a kernel's primary input data. The
  /// primary outputs get written to `flux_arr`. All 3 arrays share the shape:
  /// `(LUT::num_entries, LENz, LENy, LENx)`
  ///
  /// The index passed to axis 0 correspond to different quantities (a kernel's
  /// `LUT` specifies the precise mapping between quantities & index values).
  /// There are 2 relevant points to be mindful of:
  /// 1. for these arrays, when accessing components of vector quantities
  ///    (e.g. `LUT::velocity_i`, `LUT::velocity_j`, or `LUT::bfield_k`),
  ///    the i, j, & k components ALWAYS map to the x, y, and z components.
  /// 2. For performance purposes, the length of the arrays along axis 0 is
  ///    allowed to exceed `LUT::num_entries`. But kernels should never
  ///    access these indices (they don't map to quantities in LUT)
  /// For completeness, we provide the following examples:
  /// - `flux_arr(LUT::density,...)` & `flux_arr(LUT::velocity_k,...)`
  ///   are where the computed density & momentum_z fluxes should be stored
  /// - `prim_arr_l(LUT::density,...)` & `prim_arr_l(LUT::velocity_i,...)` hold
  ///   the reconstructed density & velocity_x values on the left side of a
  ///   cell interface
  /**@{*/
  /// array to store the computed fluxes
  const CelloView<enzo_float,4> flux_arr;
  /// array of primitives reconstructed on the left side of cell interfaces
  const CelloView<const enzo_float,4> prim_arr_l;
  /// array of primitives reconstructed on the right side of cell interfaces
  const CelloView<const enzo_float,4> prim_arr_r;
  /**@}*/

  /// @name DualEnergyArrays
  /// These arrays store outputs used in the dual-energy formalism and have the
  /// shape `(LENz, LENy, LENx)`.
  ///
  /// For any Dual-Energy compatible EOS, kernels should ALWAYS fill these
  /// arrays with the relevant values. If the dual energy formalism isn't used
  /// outside of the Riemann Solver, these are initialized with scratch-space.
  /// This is done to avoid unnecessary branching and code generation.
  /**@{*/
  /// array to store the computed flux for the specific internal energy
  ///
  /// @note
  /// When `flux_arr.shape(0) > LUT::num_entries`, this can technically alias
  /// one of `flux_arr`'s subarrays without a corresponding `LUT` entry. In
  /// practice, kernels will NEVER be affected by this.
  const CelloView<enzo_float,3> internal_energy_flux_arr;
  /// array to store the velocity (along `dim`) computed at the cell interface
  const CelloView<enzo_float,3> velocity_i_bar_arr;
  /**@}*/
};
// SPHINX-SNIPPET-KERNELCONFIG-END-INCLUDE

//----------------------------------------------------------------------

template <class KernelFunctor>
class EnzoRiemannImpl : public EnzoRiemann
{
  /// @class    EnzoRiemannImpl
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Provides implementation of approximate Riemann
  ///           Solvers
  ///
  /// @tparam KernelFunctor The functor used to specialize `EnzoRiemannImpl`.
  ///     The functor must provide a public member type called `LUT`, which is
  ///     a specialization of `EnzoRiemannLUT<InputLUT>`. It must also support
  ///     initialization from KernelConfig and provide a public `operator()`
  ///     method which accepts 3 integer indices.
  ///
  /// EnzoRiemannImpl factors out the repeated code between different
  /// approximate Riemann Solvers (e.g. HLLE, HLLC, HLLD and possibly LLF &
  /// Roe solvers).
  ///
  /// To define a new RiemannSolver using `EnzoRiemann`:
  ///   1. Define a new `KernelFunctor` (e.g. `HLLKernel`).
  ///   2. It might be useful to define an alias name for the specialization of
  ///      `EnzoRiemannImpl` that uses the new `KernelFunctor`
  ///      (e.g. `using EnzoRiemannHLLD = EnzoRiemannImpl<HLLKernel>;`).
  ///   3. Update `EnzoRiemann::construct_riemann` to construct the Riemann
  ///      Solver when the correct name is specified.
  ///   4. Update the documentation with the name of the newly available
  ///      RiemannSolver

  using LUT = typename KernelFunctor::LUT;
  using EOSStructT = typename KernelFunctor::EOSStructT;

  // Check whether KernelFunctor's operator() method has the expected signature
  // and raise an error message if it doesn't:
  static_assert(std::is_convertible<KernelFunctor&&,
                                    std::function<void(int, int, int)>>::value,
                "KernelFunctor must define the method: "
                "void KernelFunctor::operator() (int,int,int)");



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
	      const int stale_depth,
	      const str_vec_t &passive_list,
              const CelloView<enzo_float,3> * const interface_velocity) const;

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

template <class KernelFunctor>
EnzoRiemannImpl<KernelFunctor>::EnzoRiemannImpl
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

template <class KernelFunctor>
void EnzoRiemannImpl<KernelFunctor>::solve
(const EnzoEFltArrayMap &prim_map_l, const EnzoEFltArrayMap &prim_map_r,
 EnzoEFltArrayMap &flux_map, const int dim,
 const int stale_depth, const str_vec_t &passive_list,
 const CelloView<enzo_float,3> * const interface_velocity) const
{

  // Currently just going to assume that we have an Ideal Equation of State
  // (should probably test that...)
  const EnzoEOSIdeal eos_struct
    = enzo::fluid_props()->eos_variant().get<EnzoEOSIdeal>();

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
                                           eos_struct,
                                           flux_map.get_backing_array(),
                                           prim_map_l.get_backing_array(),
                                           prim_map_r.get_backing_array(),
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

template <class KernelFunctor>
void EnzoRiemannImpl<KernelFunctor>::solve_
(const KernelConfig<EOSStructT> config, const int stale_depth)
  noexcept
{
  const KernelFunctor kernel{config};

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
