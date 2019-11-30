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
  /// The source term in question for the internal energy density is formally
  /// defined as follows. We define the jth component of the interface velocity
  /// as vbar and grid width along this dimension as dx. For all cell
  /// locations (i,j,k) there is a source term, over a timestep dt, of:
  ///     -1 * dt* pressure_j * (vbar_{j+1/2} - vbar_{j+1/2}) / (a * dx_j)
  /// a corresponding term is added for each dimension (although this
  /// implementation only handles a single dimension at a time).
  ///
  /// As in Enzo's PPM method, the cell-centerred pressure is computed from the
  /// internal energy (although the internal energy and total energy should
  /// have been synchronized at the start of the timestep)
  ///
  /// The current implementation ignores the scale factor.

public:

  /// Computes the internal energy density source term contributed by the
  /// derivative of the velocity along dimension `dim` and add it to the field
  /// tracking the total change in internal_energy density.
  ///
  /// @param block holds data to be processed
  /// @param dt The time time-step overwhich to apply the source term
  /// @param prim_group holds field names of the cell-centered primitives from
  ///     the start of the timestep (but after the synchronization of the
  ///     internal energy with the total energy).
  /// @param dUcons_group Contains fields that track the accumulated change in
  ///     the integrable and passively advected quantites. The calculated
  ///     internal energy density source term will simply be added to the field
  ///     store in the "internal_energy" group. The source term for internal
  ///     energy density is computed instead of specific internal energy
  ///     because this grouping tracks changes to the conserved form of
  ///     integrable/passive quantities.
  /// @param interface_velocity_name indicates the name of field where the
  ///     value of the interface velocity along dimension `dim` is stored. This
  ///     should have been calculated by the Riemann Solver. The underlying
  ///     field should formally be defined as cell-centered to allow for it's
  ///     reuse in representing face-centered fields (excluding the exterior
  ///     faces of the block) along different dimensions
  /// @param dim Dimension along which to compute the contribution of the
  ///     internal Energy source term. Values of 0, 1, and 2 correspond to the
  ///     x, y, and z directions, respectively.
  /// @param eos Instance of the fluid's EnzoEquationOfState object
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  ///     The update using the delayed staling rate should be applied at some
  ///     time after this function call
  void calculate_source(Block *block, double dt,
			Grouping &prim_group, Grouping &dUcons_group,
			std::string interface_velocity_name, int dim,
			EnzoEquationOfState *eos, int stale_depth)
    const throw();
};

#endif /* ENZO_ENZO_SOURCE_INTERNAL_ENERGY_HPP */
