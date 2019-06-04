// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoEOSIdeal.hpp
/// @author   Matthew Abruzzo (matthewabruzzo@gmail.com)
/// @date     Thurs May 2 2019
/// @brief    [\ref Enzo] Implementation of the equation of state for an ideal
/// adiabatic gas

#ifndef ENZO_ENZO_EOS_IDEAL_HPP
#define ENZO_ENZO_EOS_IDEAL_HPP

class EnzoEOSIdeal : public EnzoEquationOfState
{

  /// @class    EnzoEOSIdeal
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates equation of state for ideal gas
  
public: // interface
  
  /// Create a new EnzoEOSIdeal object
  EnzoEOSIdeal(double gamma, double density_floor,
	       double pressure_floor) throw()
    : EnzoEquationOfState(),
      gamma_(gamma),
      density_floor_(density_floor),
      pressure_floor_(pressure_floor)
  { }

  /// Delete EnzoEOSIdeal object
  ~EnzoEOSIdeal()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoEOSIdeal);

  /// CHARM++ migration constructor for PUP::able
  EnzoEOSIdeal (CkMigrateMessage *m)
    : EnzoEquationOfState(m),
      gamma_(0.),
      density_floor_(0.),
      pressure_floor_(0.)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Computes thermal pressure
  void compute_pressure(Block *block, Grouping &cons_group,
			Grouping &prim_group, int stale_depth);

  /// Converts the cell-centered conservative quantities to primitive quantites
  void primitive_from_conservative(Block *block, Grouping &cons_group,
				   Grouping &prim_group, int stale_depth);

  /// Converts the cell-centered primitive quantities to conservative quantites
  void conservative_from_primitive(Block *block, Grouping &prim_group,
				   Grouping &cons_group, int stale_depth,
				   int reconstructed_axis);

  /// returns the density floor
  enzo_float get_density_floor() { return density_floor_; }

  /// returns the thermal pressure floor
  enzo_float get_pressure_floor() { return pressure_floor_; }

  /// apply the pressure floor to total_energy field
  void apply_floor_to_energy(Block *block, Grouping &cons_group,
			     int stale_depth);

  /// returns whether the EOS is barotropic
  bool is_barotropic() { return false; }

  /// returns adiabatic index
  enzo_float get_gamma() { return gamma_;}

  /// returns isothermal sound speed - it shouldn't be used since an ideal gas
  /// is not barotropic
  enzo_float get_isothermal_sound_speed() { return 0;}
  

private:
  /// Helper function to retrieve a field array when it is possible that a
  /// field stores reconstructed data
  EFlt3DArray retrieve_field_(EnzoFieldArrayFactory &array_factory,
			      Grouping &group, std::string group_name,
			      int index, int reconstructed_axis);
  
  /// Copies entries of the passively advected fields included by origin_group
  /// to the corresponding entries of the fields included in destination_group
  /// reconstructed_axis = -1 means that internal field shape data can be
  /// truested. Values of 0, 1, or 2 mean that the field stores reconstructed
  /// values along the x, y, or z axis and that it is actually face-centered
  void copy_passively_advected_fields_(EnzoFieldArrayFactory &array_factory,
				       Grouping &origin_group,
				       Grouping &destination_group,
				       int reconstructed_axis = -1);

protected: // attributes
  enzo_float gamma_; // adiabatic index
  enzo_float density_floor_;
  enzo_float pressure_floor_;
};

#endif /* ENZO_ENZO_EOS_IDEAL_HPP */
