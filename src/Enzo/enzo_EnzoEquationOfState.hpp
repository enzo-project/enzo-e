// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEquationOfState.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of Enzo's Enzo Equation Of State
/// abstract base class. It will be subclassed to implement specific types of
/// equations of state.

//#define RAISE_FLOOR_ERROR

#ifndef ENZO_ENZO_EQUATIONOFSTATE_HPP
#define ENZO_ENZO_EQUATIONOFSTATE_HPP
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
  /// @param reconstrable_group holds field names of reconstrable primitives
  ///     where the converted values will be stored. There is expected to be
  ///     significant overlap with the fields stored in integrable_group
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  ///
  /// For a barotropic EOS, this nominally does nothing
  /// For a non-barotropic EOS, this computes specific internal energy
  /// (unless its being tracked for the dual-energy formalism) 
  virtual void reconstructable_from_integrable(Block *block,
					       Grouping &integrable_group,
					       Grouping &reconstructable_group,
					       int stale_depth)=0;

  /// @overload
  ///
  /// Provides stale_depth with the default value of 0
  void reconstructable_from_integrable(Block *block,
				       Grouping &integrable_group,
				       Grouping &reconstructable_group)
  {
    reconstructable_from_integrable(block, integrable_group,
				    reconstructable_group, 0);
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
  /// specific internal energy
  virtual void integrable_from_reconstructable(Block *block,
					       Grouping &reconstructable_group,
					       Grouping &integrable_group,
					       int stale_depth,
					       int reconstructed_axis)=0;

  /// @overload
  ///
  /// Provides stale_depth with the default value of 0 and assumes that the
  /// fields are cell-centered
  void integrable_from_reconstructable(Block *block,
				       Grouping &reconstructable_group,
				       Grouping &integrable_group)
  { conservative_from_primitive(block, reconstr_group, integr_group, 0, -1); }

  /// Computes thermal pressure from integrable quantities
  /// 
  /// @param block holds data to be processed
  /// @param integrable_group holds field names of integrable primitives to be
  ///     used to compute thermal pressure
  /// @param pressure_name field name where the computed pressure will be
  ///     stored
  /// @param passive_scalars_group holds field names of specific passive
  ///     scalars to be (possibly used) in the calculation. This can be the
  ///     as integrable_group. These will only be used if Grackle is in use
  /// @param specific_passive_scalars indicates whether the passive scalars are
  ///     have been converted to be specific quantities
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  ///
  /// This nominally wraps EnzoComputePressure. Currently, it is primarily used
  /// to compute the pressure to determine the timestep for an integrator
  /// (before the passively advected scalars are converted from conservative
  /// quantities to primitive quantites)
  virtual void pressure_from_integrable(Block *block,
					Grouping &integrable_group,
					std::string pressure_name,
					Grouping &passive_scalars_group,
					bool specific_passive_scalars,
					int stale_depth)=0;

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
  /// This should nominally wrap EnzoComputePressure (as of now it doesn't
  /// because there is no support for MHD fields). In doing so, it should also
  /// provide grackle support
  virtual void pressure_from_reconstructable(Block *block,
					     Grouping &reconstructable_group,
					     std::string pressure_name,
					     int stale_depth,
					     int reconstructed_axis)=0;

  /// returns the density floor
  virtual enzo_float get_density_floor()=0;

  /// returns the thermal pressure floor
  virtual enzo_float get_pressure_floor()=0;

  /// apply the pressure floor to the specific total energy field. If the
  /// equation of state is barotropic, then this does nothing
  ///
  /// @param block holds data to be processed
  /// @param integrable_group holds field names of integrable primitives to be
  ///     that will be used to apply the floor (also contains the field upon
  ///     which the floor will be applied)
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  virtual void apply_floor_to_total_energy(Block *block,
					   Grouping &integrable_group,
					   int stale_depth)=0;

  /// apply the pressure floor to the specific internal energy field. If the
  /// equation of state is barotropic, then this does nothing
  ///
  /// @param block holds data to be processed
  /// @param reconstructable_group holds field names of integrable primitives
  ///     to be that will be used to apply the floor (also contains the field
  ///     upon which the floor will be applied)
  /// @param stale_depth indicates the number of field entries from the
  ///     outermost field value that the region including "stale" values (need
  ///     to be refreshed) extends over (0 means there are no "stale" values).
  ///
  /// @note If the macro RAISE_FLOOR_ERROR is defined, and a floor NEEDS to be
  /// applied, then the this function will raise an error
  virtual void apply_floor_to_internal_energy(Block *block,
					      Grouping &reconstructable_group,
					      int stale_depth)=0;

  /// returns whether the equation of state is barotropic
  virtual bool is_barotropic() = 0;

  /// returns adiabatic index - only needs to be a reasonable number of non-
  /// barotropic
  virtual enzo_float get_gamma() = 0;

  /// returns isothermal sound speed - only needs to be reasonable for a
  /// barotropic EOS
  virtual enzo_float get_isothermal_sound_speed() = 0;

  /// returns true if the dual energy formalism is being used
  virtual bool uses_dual_energy_formalism() = 0;

};

#endif /* ENZO_ENZO_EQUATIONOFSTATE_HPP */
