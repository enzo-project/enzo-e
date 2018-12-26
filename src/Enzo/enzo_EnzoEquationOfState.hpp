// This is the base class defining the interface for objects encapsulating the
// equation of state of the fluid. It will be subclassed to represent different
// equations of state (e.g. ideal gas, polytropic EOS, isothermal, etc.)
//
// Current plan: track instances within the hydro method.
//
// Several of the methods, are expecting C-style arrays. Current plan is that
// the EquationOfState will know how many variables to expect (e.g. if using
// dual energy formalism)

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

  /// Virtual destructor
  virtual ~EnzoEquationOfState()
  { }

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

  // Converts the cell-centered conservative quantities to primitive quantites
  virtual void primitive_from_conservative(Block * block, Grouping &cons_group,
  					   Grouping &prim_group)=0;

  virtual void conservative_from_primitive(enzo_float *prim, enzo_float *cons);

  // Computes magnetic pressure from primitives
  virtual enzo_float mag_pressure_from_primitive(enzo_float *prim_vals)=0;
  
  // Computes the thermal sound speed
  virtual enzo_float sound_speed (enzo_float *prim_vals)=0;

  // computes the fast magnetosonic speed
  virtual enzo_float fast_magnetosonic_speed (enzo_float *prim_vals) =0;

  // returns adiabatic index
  virtual enzo_float get_gamma() = 0;

  // returns the density floor
  virtual enzo_float get_density_floor()=0;

  // returns the thermal pressure floor
  virtual enzo_float get_pressure_floor()=0;

  /// Name of this equation of state type
  virtual std::string type () const
  { return "undefined"; }

};

#endif /* ENZO_ENZO_EQUATIONOFSTATE_HPP */
