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
  static void check_floor(enzo_float floor, bool density){
    if (density){
      // if density = 0, NaNs will arise when converting momentum to velocity
      ASSERT("EnzoEquationOfState::check_floor",
	     "The density floor must be greater than 0",
	     floor > 0);
    } else {
      ASSERT("EnzoEquationOfState::check_floor",
	     "The floor must be greater than or equal to 0",
	     floor >= 0);
    }
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

  /// Computes thermal pressure - will ideally wrap EnzoComputePressure
  virtual void compute_pressure(Block *block, Grouping &cons_group,
				Grouping &prim_group, int stale_depth)=0;

  /// Converts the cell-centered conservative quantities to primitive quantites
  virtual void primitive_from_conservative(Block *block, Grouping &cons_group,
  					   Grouping &prim_group,
					   int stale_depth)=0;

  // Provides stale_depth with the default value of 0 
  void primitive_from_conservative(Block *block, Grouping &cons_group,
  				   Grouping &prim_group)
  { primitive_from_conservative(block, cons_group, prim_group, 0); }

  /// Converts the cell-centered primitive quantities to conservative quantites
  /// @param block holds data to be processed
  /// @param prim_group holds field names of primitives values to convert
  /// @param cons_group holds field names of conservatives where the converted
  /// values will be stored
  /// @param stale_depth indicates the number of field entries from the
  /// outermost field value that the region including "stale" values (need to
  /// be refreshed) extends over (0 means there are no "stale" values).
  /// @param reconstructed_axis - parameter that optionally indicates that the
  /// primitives have been reconstructed (if the same array is used multiple
  /// axes, then the fields are internally stored as cell-centered which is
  /// than an array of face-centered quantites). A value of -1 means that
  /// they are cell-centered. A value of 0, 1, or 2 means that the fields were
  /// reconstructed and they only contain valid values at x, y, or z faces
  virtual void conservative_from_primitive(Block *block, Grouping &prim_group,
  					   Grouping &cons_group,
					   int stale_depth,
					   int reconstructed_axis)=0;

  // Provides stale_depth with the default value of 0 and assumes cell-centered
  void conservative_from_primitive(Block *block, Grouping &cons_group,
  				   Grouping &prim_group)
  { conservative_from_primitive(block, cons_group, prim_group, 0, -1); }

  /// returns the density floor
  virtual enzo_float get_density_floor()=0;

  /// returns the thermal pressure floor
  virtual enzo_float get_pressure_floor()=0;

  /// apply the pressure floor to total_energy field
  virtual void apply_floor_to_energy(Block *block, Grouping &cons_group,
				     int stale_depth)=0;

  /// returns whether the equation of state is barotropic
  virtual bool is_barotropic() = 0;

  /// returns adiabatic index - only needs to be a reasonable number of non-
  /// barotropic
  virtual enzo_float get_gamma() = 0;

  /// returns isothermal sound speed - only needs to be reasonable for a
  /// barotropic EOS
  virtual enzo_float get_isothermal_sound_speed() = 0;

};

#endif /* ENZO_ENZO_EQUATIONOFSTATE_HPP */
