#ifndef ENZO_ENZO_EOS_IDEAL_HPP
#define ENZO_ENZO_EOS_IDEAL_HPP

// This should almost certainly wrap an instance of ComputePressure (and
// presumably Compute Temperature)
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

  // Computes thermal pressure
  void compute_pressure(Block *block, Grouping &cons_group,
			Grouping &prim_group);

  // Converts the cell-centered conservative quantities to primitive quantites
  void primitive_from_conservative(Block *block, Grouping &cons_group,
				   Grouping &prim_group);

  void conservative_from_primitive(Block *block, Grouping &prim_group,
				   Grouping &cons_group);

  // Computes magnetic pressure from primitives
  enzo_float mag_pressure_from_primitive(flt_map &prim_vals){
    return 0.5 * (prim_vals["bfield_i"]*prim_vals["bfield_i"]+
		  prim_vals["bfield_j"]*prim_vals["bfield_j"]+
		  prim_vals["bfield_k"]*prim_vals["bfield_k"]);}

  // Computes the thermal sound speed
  enzo_float sound_speed(flt_map &prim_vals){
    return std::sqrt(gamma_*prim_vals["pressure"]/prim_vals["density"]); }

  // computes the fast magnetosonic speed along dimension i
  enzo_float fast_magnetosonic_speed(flt_map &prim_vals);

  // returns adiabatic index
  enzo_float get_gamma(){
    return gamma_;}

  // returns the density floor
  enzo_float get_density_floor(){
    return density_floor_;}

  // returns the thermal pressure floor
  enzo_float get_pressure_floor(){
    return pressure_floor_;}

  // apply the pressure floor to total_energy field
  void apply_floor_to_energy(Block *block, Grouping &cons_group);

protected: // attributes
  enzo_float gamma_; // adiabatic index
  enzo_float density_floor_;
  enzo_float pressure_floor_;
};

#endif /* ENZO_ENZO_EOS_IDEAL_HPP */
