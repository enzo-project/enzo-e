// This is the base class defining the interface for objects encapsulating the
// equation of state of the fluid. It will be subclassed to represent different
// equations of state (e.g. ideal gas, polytropic EOS, isothermal, etc.)
//
// Current plan: track instances within the hydro method.

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

  /// CHARM++ PUP::able declaration
  //PUPable_abstract(EnzoEquationOfState);

  /// CHARM++ migration constructor for PUP::able
  //EnzoEquationOfState (CkMigrateMessage *m)
  //  : PUP::able(m)
  //{  }

  /// Virtual destructor
  virtual ~EnzoEquationOfState()
  {  }

  /// CHARM++ Pack / Unpack function
  //void pup (PUP::er &p)
  //{
  //  PUP::able::pup(p);
  //}

  // Converts the cell-centered conservative quantities to primitive quantites
  virtual void conservative_from_primitive (Block * block,
					    std::vector<int> &cons_ids,
  					    std::vector<int> &prim_ids)=0;
	    
  // Converts primative quantities to conservative quantites
  virtual void primitive_from_conservative (Block * block,
					    std::vector<int> &cons_ids,
  					    std::vector<int> &prim_ids)=0;

  // Computes the thermal sound speed
  virtual enzo_float sound_speed ()=0;

  // computes the fast magnetosonic speed
  virtual enzo_float fast_magnetosonic_speed () =0;

  /// Name of this equation of state type
  virtual std::string type () const
  { return "undefined"; }

};

#endif /* ENZO_ENZO_EQUATIONOFSTATE_HPP */
