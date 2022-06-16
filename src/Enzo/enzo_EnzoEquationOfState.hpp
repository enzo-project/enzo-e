// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEquationOfState.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Equation Of State
/// abstract base class. It will be subclassed to implement specific types of
/// equations of state.

//#define RAISE_FLOOR_ERROR

#ifndef ENZO_ENZO_EQUATIONOFSTATE_HPP
#define ENZO_ENZO_EQUATIONOFSTATE_HPP

// Among it's EOS-related responsbilities, EnzoEquationOfState, is responsible
// for the application of the Dual Energy Formalism, (when specified for
// non-barotropic equations of state). Currently, implementations of the Dual
// Energy Formalism are expected to more closely resemble implementations in
// Enzo's Runge-Kutta and MHD with Constrained Transport integrators. These
// exhibit 3 main differences from the original conception (implemented in
// Enzo's ppm integrator):
//     1. internal energy is always used to compute pressure. In the original
//        conception, pressure could be computed from total energy or
//        internal energy (the decision was independent of synchronization).
//     2. Unlike the original conception, both pressure and internal energy are
//        not reconstructed separately. Implementations are currently expected
//        to just reconstruct pressure and compute internal energy from the
//        reconstructed quantities.
//     3. Synchronization of the total and internal energies is a local
//        operation that doesn't require knowledge cell neighbors. In the
//        original conception, knowledge of the immediate neighbors had been
//        required (thus, each synchronization incremented the stale depth).
//
// To allow synchronization to data from neigboring cells, an additional method
// would be required that indicates the staling_rate of the synchronization.

class EnzoEquationOfState : public PUP::able
{
  /// @class    EnzoEquationOfState
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates equation of state of fluid

public: // interface

  /// Create a new EnzoEquationOfState
  EnzoEquationOfState() throw()
  {}

  ~EnzoEquationOfState()
  { }

  /// Checks the validity of floor values for the EquationOfState
  static void check_floor(enzo_float floor){
    // if density = 0, NaNs will arise when converting momentum to velocity
    // if pressure = 0, then sound speed will be equal to 0 (possibly causing
    // time-step calculation problems)
    ASSERT("EnzoEquationOfState::check_floor",
	   "The density and pressure floors must be greater than 0",
	   floor > 0);
  }

  /// Applies primitive floor. This function has been factored out to allow for
  /// more easily debugging cases when the floor is unnecesarily applied.
  inline static enzo_float apply_floor(const enzo_float value,
                                       const enzo_float floor){
    enzo_float out;
    #ifdef RAISE_FLOOR_ERROR
    ASSERT("EnzoEquationOfState", "Should not need to apply primitive floor.",
	   value >= floor);
    out = value;
    #else
    out = std::max(value,floor);
    #endif
    return out;
  }

  /// CHARM++ PUP::able declaration
  PUPable_abstract(EnzoEquationOfState);

  /// CHARM++ migration constructor for PUP::able
  EnzoEquationOfState (CkMigrateMessage *m)
    : PUP::able(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
  }

  /// Converts integration quantities to primitives
  ///
  /// @param[in]  integration_map Map holding integration quantities that are
  ///     to be converted. Passive scalars in this map are expected to be in
  ///     conserved form (they are densities).
  /// @param[out] primitive_map Map holding arrays where the computed
  ///     primitive data is to be stored. Passive scalars in this map are
  ///     expected to be in specific form (they are mass fractions).
  /// @param[in]  stale_depth indicates the current stale_depth for the
  ///     supplied cell-centered quantities
  /// @param[in]  passive_list A list of keys for passive scalars. These keys
  ///     will be used to determine which quantities will be copied from the
  ///     integration_map to the primitive_map.
  ///
  /// Non-passive scalar quantities appearing in both `integration_map` and
  /// `primitive_map` are simply deepcopied and passive scalar quantities are
  /// converted from conserved-form to specific form. For a non-barotropic EOS,
  /// this also computes pressure.
  ///
  /// @note
  /// It's a not obvious to me that the EOS should necessarily be responsible
  /// for this operation.
  virtual void primitive_from_integration
  (const EnzoEFltArrayMap &integration_map, EnzoEFltArrayMap &primitive_map,
   const int stale_depth, const str_vec_t &passive_list) const =0;

  /// Computes thermal pressure from integration quantities
  /// 
  /// @param[in]  integration_map Map holding integration quantities that are
  ///     used to compute the pressure. This should include all necessary
  ///     passively advected quantities in conserved form.
  /// @param[out] pressure Array where the thermal pressure is to be stored
  /// @param[in]  stale_depth indicates the current stale_depth for the
  ///     supplied cell-centered quantities
  ///
  /// This nominally should wrap EnzoComputePressure. But at the time of
  /// writing, it doesn't actually wrap EnzoComputePressure
  virtual void pressure_from_integration
  (const EnzoEFltArrayMap &integration_map,
   const CelloArray<enzo_float, 3> &pressure,
   const int stale_depth) const = 0;

  /// returns the density floor
  virtual enzo_float get_density_floor() const = 0;

  /// returns the thermal pressure floor
  virtual enzo_float get_pressure_floor() const = 0;

  /// applies the pressure floor to the specific total energy field. If using
  /// the dual-energy formalism, it is also applied to the internal energy
  /// and it synchronize the internal energy and total energy fields. If the
  /// EOS is barotropic, this does nothing.
  ///
  /// @param[in,out] integration_map Map holding integration quantities that
  ///     will be used to apply the floor. It must also include a
  ///     "total_energy" entry (unless the EOS is barotropic) upon which the
  ///     floor is applied.
  /// @param[in]     stale_depth indicates the current stale_depth for the
  ///     supplied cell-centered quantities
  ///
  /// Unlike the initial conception of the dual-energy formalism (or the
  /// version used in Enzo's ppm integrator) this assumes that synchronization
  /// is a local operation that doesn't require data about neighboring cells
  /// (similar to the implementation of the dual energy formalsim in Enzo's
  /// Runge Kutta and MHD with Constrained Transport solvers).
  virtual void apply_floor_to_energy_and_sync(EnzoEFltArrayMap &integration_map,
                                              const int stale_depth) const = 0;

  /// returns whether the equation of state is barotropic
  virtual bool is_barotropic() const = 0;

  /// returns adiabatic index - only needs to be a reasonable number of non-
  /// barotropic
  virtual enzo_float get_gamma() const = 0;

  /// returns isothermal sound speed - only needs to be reasonable for a
  /// barotropic EOS
  virtual enzo_float get_isothermal_sound_speed() const = 0;

  /// returns true if the dual energy formalism is being used
  virtual bool uses_dual_energy_formalism() const = 0;

};

#endif /* ENZO_ENZO_EQUATIONOFSTATE_HPP */
