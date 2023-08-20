// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSourceInternalEnergy.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Mon December 1 2019
/// @brief    [\ref Enzo] Declaration of Enzo's SourceInternalEnergy class.
/// This computes the source term in the internal energy due to the velocity
/// gradient along a given dimension.

// This is computes the source terms in a dimensionally split manner so it can
// be used for dimensionally split methods (like PPM) in the future. For this
// reason and because it should be applied in predictor-corrector methods (like
// VL+CT) at partial timesteps, this has not been operator split off into a
// separate method
//
// Like in Enzo's PPM method, this computes spatial derivatives in velocity
// using the interface velocities (as opposed to Enzo's Runge-Kutta method
// which uses cell-centered velocities).
//
// The calculation of the source term is not performed within the approximate
// Riemann Solver as it requires knowledge of the cell-centered
//
// As additional source (or loss) terms are introduced, it may make sense to
// introduce a class hierarchy to simply the implementation/application of
// source terms

#ifndef ENZO_ENZO_SOURCE_INTERNAL_ENERGY_HPP
#define ENZO_ENZO_SOURCE_INTERNAL_ENERGY_HPP

class EnzoSourceInternalEnergy
{
  /// @class    EnzoSourceInternalEnergy
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the calculation of the internal energy
  ///           source term along a given dimension.
  ///
  /// For convenience, we split the encapsulated source term up by into 3
  /// pieces; one for each dimension (this class needs to be called separately
  /// for each dimension). We define the portion of the source term for the ith
  /// dimension that get's applied at every cell locations as:
  /// @code{.unparsed}
  ///    -1 * dt * pressure[i] * (vbar[i+1/2] - vbar[i-1/2]) / (a * cellwidth_i)
  /// @endcode
  /// In the above equation:
  ///   * the bracketted variables indicate the relative positions of variables
  ///     on the grid along the ith dimension (values of `i` map to a
  ///     cell-centered position). You can assume that the omitted indices
  ///     along the other dimensions are always cell-centered.
  ///   * `vbar` is the ith component of the velocity at the cell interface
  ///   * `a` is the current scale-factor.
  ///   * `cellwidth_i` is the comoving cell width along the ith dimension. In
  ///     practice, this calculation always expects to be passed the proper
  ///     cell-width, which is always equal to `a * cellwidth_i`.
  ///
  /// This source-term is compatible with both non-cosmological and
  /// cosmological simulations. In cosmological calculations the internal
  /// energy density has another source term, but that get's handled separately
  /// by `EnzoMethodComovingExpansion`.
  ///
  /// As in Enzo's PPM method, the cell-centered pressure is computed from the
  /// internal energy (although the internal energy and total energy should
  /// have been synchronized at the start of the timestep)

public:

  /// Computes the internal energy density source term contributed by the
  /// derivative of the velocity along dimension `dim` and add it to the array
  /// tracking the total change in internal_energy density.
  ///
  /// @param[in]  dim Dimension along which to compute the contribution of the
  ///     internal Energy source term. Values of 0, 1, and 2 correspond to the
  ///     the x, y, and z directions, respectively.
  /// @param[in]  dt The time time-step overwhich to apply the source term
  /// @param[in]  proper_cell_width The cell width along dimension `dim`. This
  ///     should always be a proper distance (even in cosmological sims).
  /// @param[in]  prim_map Map holding the values of the cell-centered
  ///     primitives from the start of the (partial) timestep (but after the
  ///     synchronization of the internal energy with the total energy).
  ///     The cell-centered pressure is loaded from this.
  /// @param[out] dUcons_map Map of arrays where the net changes to the
  ///     integration quantities (including passively advected scalars) are
  ///     accumulated. The internal energy density source term is simply added
  ///     to the array associated with the "internal_energy" key. The source
  ///     term for internal energy density is computed instead of specific
  ///     internal energy because this map accumulates the net change to the
  ///     conserved form of the integration/passive quantities.
  /// @param[in]  interface_velocity Array storing the values of the interface
  ///     velocity along dimension `dim` is stored (it should be computed by
  ///     the Riemann Solver). This should have the same dimensions as the
  ///     arrays stored in prim_map/dUcons_map, except along dimension `dim`.
  ///     Along that dimension, this array should hold one fewer value.
  /// @param[in]  stale_depth indicates the current stale_depth for the
  ///     supplied cell-centered quantities. The update using the
  ///     reconstructor's delayed_staling_rate should be applied at some
  ///     time after this function call.
  void calculate_source(const int dim, const double dt,
			const enzo_float proper_cell_width,
                        const EnzoEFltArrayMap &prim_map,
                        EnzoEFltArrayMap &dUcons_map,
                        const CelloView<const enzo_float,3> &interface_velocity,
			const int stale_depth)
    const throw();
};

#endif /* ENZO_ENZO_SOURCE_INTERNAL_ENERGY_HPP */
