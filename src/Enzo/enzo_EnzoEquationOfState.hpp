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
  static enzo_float apply_floor(const enzo_float value, const enzo_float floor){
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


  /// Converts cell-centered integrable primitives to reconstructable primitives
  ///
  /// @param block holds data to be processed
  /// @param integrable_group holds field names of integrable primitive values
  ///     to convert
  /// @param reconstructable_group holds field names of reconstructable
  ///     primitives where the converted values will be stored. There is
  ///     expected to be significant overlap with the fields stored in
  ///     integrable_group
  /// @param conserved_passive_group contains the names of the fields holding
  ///     the passively advected scalars in conserved form (note that while the
  ///     integrable grouping may also contain groups of passive scalar fields,
  ///     those fields hold the passive scalars in specific form - which are
  ///     never used in this calculation). These are provided for Grackle's use
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  ///
  /// For a barotropic EOS, this nominally does nothing
  /// For a non-barotropic EOS, this computes pressure 
  virtual void reconstructable_from_integrable
  (Block *block, Grouping &integrable_group, Grouping &reconstructable_group,
   Grouping &conserved_passive_group, int stale_depth) const = 0;

  /// @overload
  ///
  /// Provides stale_depth with the default value of 0
  void reconstructable_from_integrable(Block *block,
				       Grouping &integrable_group,
				       Grouping &reconstructable_group,
				       Grouping &conserved_passive_group) const
  {
    reconstructable_from_integrable(block, integrable_group,
				    reconstructable_group,
				    conserved_passive_group, 0);
  }

  /// Converts reconstructable primitives to integrable primitives
  ///
  /// @param block holds data to be processed
  /// @param reconstructable_group holds field names of reconstructable
  ///     primitive values to convert
  /// @param integrable_group holds field names of integrable primitives where
  ///     the converted values will be stored. There is expected to be
  ///     significant overlap with the fields stored in reconstructable_group
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  /// @param reconstructed_axis - parameter that optionally indicates that the
  ///     reconstructable primitives have been reconstructed at face-centers (if
  ///     the same grouping is used for multiple axes, then the fields are
  ///     internally stored as cell-centered rather than an array of
  ///     face-centered quantites). A value of -1 means that they are
  ///     cell-centered. A value of 0, 1, or 2 means that the fields were
  ///     reconstructed and they only contain valid values at x, y, or z faces
  ///
  /// For a barotropic EOS, this nominally does nothing
  /// For a non-barotropic EOS, this computes specific total energy from
  /// pressure. If using the dual energy formalism, it will also compute the
  /// internal energy from the pressure
  virtual void integrable_from_reconstructable(Block *block,
					       Grouping &reconstructable_group,
					       Grouping &integrable_group,
					       int stale_depth,
					       int reconstructed_axis) const =0;

  /// @overload
  ///
  /// Provides stale_depth with the default value of 0 and assumes that the
  /// fields are cell-centered
  void integrable_from_reconstructable(Block *block,
				       Grouping &reconstructable_group,
				       Grouping &integrable_group) const
  { integrable_from_reconstructable(block, reconstructable_group,
				    integrable_group, 0, -1); }

  /// Computes thermal pressure from integrable quantities
  /// 
  /// @param block holds data to be processed
  /// @param integrable_group holds field names of integrable primitives to be
  ///     used to compute thermal pressure
  /// @param pressure_name field name where the computed pressure will be
  ///     stored
  /// @param conserved_passive_group contains the names of the fields holding
  ///     the passively advected scalars in conserved form (note that while the
  ///     integrable grouping may also contain groups of passive scalar fields,
  ///     those fields hold the passive scalars in specific form - which are
  ///     never used in this calculation). These are provided for Grackle's use
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  ///
  /// This nominally should wrap EnzoComputePressure. At the time of writing,
  /// (Grackle not yet supported), it doesn't actually wrap EnzoComputePressure
  virtual void pressure_from_integrable(Block *block,
					Grouping &integrable_group,
					std::string pressure_name,
					Grouping &conserved_passive_group,
					int stale_depth) const =0;

  /// Computes thermal pressure from reconstructable quantities (nominally
  /// after reconstruction)
  ///
  /// @param block holds data to be processed
  /// @param reconstructable_group holds field names of reconstructable
  ///     primitives to be used to compute thermal pressure
  /// @param pressure_name field name where the computed pressure will be
  ///     stored
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  /// @param reconstructed_axis - parameter that optionally indicates that the
  ///     reconstructable primitives have been reconstructed at face-centers (if
  ///     the same grouping is used for multiple axes, then the fields are
  ///     internally stored as cell-centered rather than an array of
  ///     face-centered quantites). A value of -1 means that they are
  ///     cell-centered. A value of 0, 1, or 2 means that the fields were
  ///     reconstructed and they only contain valid values at x, y, or z faces
  ///
  /// For a non-barotropic EOS, pressure is considered a reconstructable
  /// quantity. In that case, if the pressure field in reconstructable_group
  /// matches pressure_name, nothing happens. If the field names do not match,
  /// then values are simply copied
  virtual void pressure_from_reconstructable(Block *block,
					     Grouping &reconstructable_group,
					     std::string pressure_name,
					     int stale_depth,
					     int reconstructed_axis) const = 0;

  /// returns the density floor
  virtual enzo_float get_density_floor() const = 0;

  /// returns the thermal pressure floor
  virtual enzo_float get_pressure_floor() const = 0;

  /// apply the pressure floor to the specific total energy field and (if using
  /// the dual-energy formalism) synchronize the internal energy and total
  /// energy fields. If the EOS is barotropic, this does nothing.
  ///
  /// @param block holds data to be processed
  /// @param integrable_group holds field names of integrable primitives to be
  ///     that will be used to apply the floor (also contains the field upon
  ///     which the floor will be applied)
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  ///
  /// Unlike the initial conception of the dual-energy formalism (or the
  /// version used in Enzo's ppm integrator) this assumes that synchronization
  /// is a local operation that doesn't require data about neighboring cells
  /// (similar to the implementation of the dual energy formalsim in Enzo's
  /// Runge Kutta and MHD with Constrained Transport solvers).
  virtual void apply_floor_to_energy_and_sync(Block *block,
					      Grouping &integrable_group,
					      int stale_depth) const = 0;

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
