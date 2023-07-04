// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSourceInternalEnergy.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Fri May 6 2022
/// @brief    [\ref Enzo] Declaration of Enzo's SourceGravity class. This
/// computes the momentum and energy source terms from gravity (or any other
/// contributors to the acceleration field).


#ifndef ENZO_ENZO_SOURCE_GRAVITY_HPP
#define ENZO_ENZO_SOURCE_GRAVITY_HPP

class EnzoSourceGravity
{
  /// @class    EnzoSourceGravity
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates the calculation of the gravity source
  ///           term (for hydro solvers)
  ///
  /// This is a pretty simple approach that was adopted from Enzo's Runge-Kutta
  /// solver.
  ///
  /// Currently, this assumes a 3D simulation. In the future, we may want to
  /// make this a configurable parameter that is passed to the constructor.

public:

  /// Computes the internal energy density source term contributed by the
  /// derivative of the velocity along dimension `dim` and add it to the array
  /// tracking the total change in internal_energy density.
  ///
  /// The gravity source terms are added to the arrays in `dUcons_map` with the
  /// associated with the "velocity_x", "velocity_y", "velocity_z", and
  /// "total_energy" keys. The source terms are computed for each component
  /// of the momentum density and the total energy density (rather than for the
  /// velocity components and specific total energy) because `dUcons_map`
  /// accumulates the net change to the conserved form of the
  /// integration/passive quantities.
  ///
  /// If the equation of state is barotropic, then this will not compute the
  /// source term for the "total_energy".
  ///
  /// @param[in]  dt The time time-step overwhich to apply the source term
  /// @param[in]  prim_map Map holding the cell-centered density and velocity
  ///     components from the start of the timestep.
  /// @param[out] dUcons_map Map of arrays where the net changes to the
  ///     integration quantities (including passively advected scalars) are
  ///     accumulated.
  /// @param[in]  accel_map Map of arrays that holds the acceleration
  ///     components. The expected keys are "acceleration_x", "acceleration_y",
  ///     and "acceleration_z". These arrays should have the same dimensions as
  ///     the arrays stored in prim_map/dUcons_map.
  /// @param[in]  stale_depth indicates the current stale_depth for the
  ///     supplied cell-centered quantities. The update using the
  ///     reconstructor's delayed_staling_rate should be applied at some
  ///     time after this function call.
  void calculate_source(const double dt,
                        const EnzoEFltArrayMap &prim_map,
                        EnzoEFltArrayMap &dUcons_map,
                        const EnzoEFltArrayMap &accel_map,
			const int stale_depth) const noexcept;
};

#endif /* ENZO_ENZO_SOURCE_GRAVITY_HPP */
