// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemann.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of the Riemann Solver abstract base
/// class. This class should be subclassed to implement various riemann solvers.
///
/// This header also serves as the public header of the EnzoRiemann sublibrary.
/// No other headers from this library should be included outside of this
/// directory

#ifndef ENZO_ENZO_RIEMANN_HPP
#define ENZO_ENZO_RIEMANN_HPP


// The conventional approach of having EnzoRiemann inherit from PUP::able
// (and making this library into a charm++ module) complicates the build-system
// to some extent. Instead, we have defined the global ``operator|``
// function to encapsulate the PUP functionality

class EnzoRiemann
{
  /// @class    EnzoRiemann
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate approximate Riemann Solvers

public: // interface

  struct FactoryArgs{
    /// @class    EnzoRiemannFactoryArgs
    /// @ingroup  Enzo
    /// @brief    [\ref Enzo] Stores arguments for EnzoRiemann's factory method

    /// The name of the Riemann solver to use. Valid names include
    /// "hll", "hlle", and "hlld"
    std::string solver;
    /// Indicates whether magnetic fields are present
    bool mhd;
    /// Indicates if the specific internal energy is an integration quantity
    bool internal_energy;
  };

  /// Factory method for constructing the EnzoRiemann object.
  static EnzoRiemann* construct_riemann(const FactoryArgs& factory_args)
    noexcept;

  /// friend function that implements CHARM++ Pack / Unpack functionality
  friend inline void operator|(PUP::er &p, EnzoRiemann*& riemann_ptr){
    bool is_nullptr = (riemann_ptr == nullptr);
    p|is_nullptr;
    if (!is_nullptr){
      FactoryArgs factory_args;
      if (!p.isUnpacking()){ factory_args = riemann_ptr->factory_args_; }

      // NOTE: change this function when the members of FactoryArgs change
      p | factory_args.solver;
      p | factory_args.mhd;
      p | factory_args.internal_energy;

      if (p.isUnpacking()){
        riemann_ptr = EnzoRiemann::construct_riemann(factory_args);
      }
    } else {
      riemann_ptr = nullptr;
    }
  }

  EnzoRiemann(const FactoryArgs& factory_args) noexcept
    : factory_args_(factory_args)
  {}

  /// Virtual destructor
  virtual ~EnzoRiemann()
  {}

  /// Computes the Riemann Fluxes for each conserved field along a given
  /// dimension, dim
  /// @param[in]     priml_map,primr_map Maps of arrays holding the left/right
  ///     reconstructed face-centered primitives. A list of the expected
  ///     keys (and key-order) is provided by the `primitive_quantity_keys`
  ///     method. These should be face-centered along `dim` (without having
  ///     values on the exterior faces of the block) and cell-centered along
  ///     the other dimensions.
  /// @param[out]    flux_map Holds arrays where the calculated fluxes for the
  ///     integration quantities will be stored. This is expected to have the
  ///     keys and key-order specified by the `integration_quantity_keys`
  ///     method. The arrays should be face-centered along `dim` (without
  ///     having values on the exterior faces of the block)
  /// @param[in]     dim Dimension along which to compute Riemann fluxes.
  ///     Values of 0, 1, and 2 correspond to the x, y, and z directions.
  /// @param[in]     eos Instance of the fluid's EnzoEquationOfState object
  /// @param[in]     stale_depth indicates the current stale_depth.
  /// @param[in]     passive_list A list of keys for passive scalars.
  /// @param[in,out] interface_velocity Pointer to an array to hold the
  ///     computed component of the velocity at the cell interfaces along
  ///     `dim` (the array should not include exterior faces of the block and
  ///     should be cell-centered along other dimensions). This quantity is
  ///     used to compute the internal energy source term (needed under the
  ///     dual energy formalism). If the value is `nullptr`, then the interface
  ///     velocity is not stored in the array.
  ///
  /// @note This function expects that the keys within `priml_map` and
  /// `primr_map` ordered such that the passive scalar keys occur after the
  /// keys specified by `primitive_quantity_keys()`. Likewise, in ``flux_map``,
  /// the passive scalar keys should occur after the keys specified by
  /// `integration_quantity_keys()`.
  ///
  /// @note This function expects that calling the `contiguous_arrays()`
  /// instance method for `prim_map_l`, `prim_map_r`, and `flux_map` will
  /// return `true` in each case.
  virtual void solve
  (const EnzoEFltArrayMap &prim_map_l, const EnzoEFltArrayMap &prim_map_r,
   EnzoEFltArrayMap &flux_map, const int dim, const EnzoEquationOfState *eos,
   const int stale_depth, const str_vec_t &passive_list,
   const CelloArray<enzo_float,3> * const interface_velocity) const = 0;

  /// Return the expected keys (and key-order) that the `solve` method expects
  /// the `flux_map` argument to have (i.e. these correspond to the fluxes of
  /// the actively advected quantities that the solver computes)
  virtual const std::vector<std::string> integration_quantity_keys()
    const noexcept = 0;

  /// Return the expected keys (and key-order) that the `solve` method expects
  /// the `priml_map` and `primr_map` arguments to have (i.e. these are the
  /// keys for the primitives that are required to compute the flux)
  virtual const std::vector<std::string> primitive_quantity_keys()
    const noexcept = 0;

private:

  /// debugging method (checks that the order of keys matches expectations
  void check_key_order_(const EnzoEFltArrayMap &map, bool prim,
			const str_vec_t &passive_list) const noexcept;

private:
  /// this is only stored to facilitate (de-)serialization
  FactoryArgs factory_args_;
};

#endif /* ENZO_ENZO_RIEMANN_HPP */
