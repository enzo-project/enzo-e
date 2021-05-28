// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoRiemann.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of the Riemann Solver abstract base
/// class. This class should be subclassed to implement various riemann solvers.

#ifndef ENZO_ENZO_RIEMANN_HPP
#define ENZO_ENZO_RIEMANN_HPP



class EnzoRiemann : public PUP::able
{
  /// @class    EnzoRiemann
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate approximate Riemann Solvers

public: // interface

  /// Factory method for constructing the EnzoRiemann object. (The signature
  /// may need to be modified as additional physics get added)
  ///
  /// @param solver The name of the Riemann solver to use. Valid names include
  ///     "hll", "hlle", and "hlld"
  /// @param mhd Indicates whether magnetic fields are present
  /// @param internal_energy Indicates whether the internal energy is an
  ///     integration quantity
  static EnzoRiemann* construct_riemann(const std::string& solver, const bool mhd,
                                        const bool internal_energy);

  EnzoRiemann() noexcept
  {}

  /// Virtual destructor
  virtual ~EnzoRiemann()
  {}

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoRiemann);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemann (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
  }

  /// Computes the Riemann Fluxes for each conserved field along a given
  /// dimension, dim
  /// @param[in]     priml_map,primr_map Maps of arrays holding the left/right
  ///     reconstructed face-centered primitives. A list of the expected
  ///     primitive quantities is provided by the `primitive_quantities` method.
  ///     These should be face-centered along `dim` (without having values on
  ///     the exterior faces of the block) and cell-centered along the other
  ///     dimensions.
  /// @param[out]    flux_map Holds arrays where the calculated fluxes for the
  ///     integration quantities will be stored. Fluxes are computed for each
  ///     integration quantities in the list provided by the
  ///     ``integration_quantities`` method. The arrays should be face-centered
  ///     along `dim` (without having values on the exterior faces of the block)
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
  /// @note It's alright for arrays in `priml_map` and `primr_map` to have the
  /// shapes of cell-centered arrays. In this case, the function effectively
  /// treats such arrays as if their `subarray` method were invoked, where
  /// `CSlice(0,-1)` is specified for the `dim` axis and
  /// `CSlice(nullptr,nullptr)` is specified for other axes. This logic also
  /// applies to the other arrays passed as arguments.
  virtual void solve
  (EnzoEFltArrayMap &prim_map_l, EnzoEFltArrayMap &prim_map_r,
   EnzoEFltArrayMap &flux_map, const int dim, const EnzoEquationOfState *eos,
   const int stale_depth, const str_vec_t &passive_list,
   const EFlt3DArray *interface_velocity) const = 0;

  /// Return the names of the actively advected integration quantities for
  /// which fluxes will be computed.
  virtual const std::vector<std::string> integration_quantities()
    const noexcept = 0;

  /// Return the names of the (non passive scalar) primitives that are required
  /// to compute the flux.
  virtual const std::vector<std::string> primitive_quantities()
    const noexcept = 0;
};

#endif /* ENZO_ENZO_RIEMANN_HPP */
