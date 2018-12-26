#ifndef ENZO_ENZO_RIEMANN_HLLE_HPP
#define ENZO_ENZO_RIEMANN_HLLE_HPP
class EnzoRiemannHLLE : public EnzoRiemann
{

  /// @class    EnzoRiemannHLLE
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulates HLLE approximate Riemann Solver

public: // interface
  /// Create a new EnzoRiemannHLLE object
  EnzoRiemannHLLE() throw()
  {}

  /// Virtual destructor
  virtual ~EnzoRiemannHLLE()
  {  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoRiemannHLLE);

  /// CHARM++ migration constructor for PUP::able
  EnzoRiemannHLLE (CkMigrateMessage *m)
    : EnzoRiemann(m)
  {  }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    PUP::able::pup(p);
  };

  /// Solve the Riemann Problem - dim (0, 1, or 2) tells the solver which
  /// dimension to compute fluxes along and indicates the dimension along which
  /// the reconstructed primitive values are face-centered
  void solve (Block *block, Grouping &priml_group, Grouping &primr_group,
	      Grouping &flux_group, int dim, EnzoEquationOfState *eos);
  
protected: // methods

  // Computes the Wave Speeds bp and bm
  // The interface is subject to change
  void wave_speeds_ (enzo_float *wl, enzo_float *wr, enzo_float *Ul,
		     enzo_float *Ur, enzo_float mag_p_l, enzo_float mag_p_r,
		     EnzoEquationOfState *eos, enzo_float *bp, enzo_float *bm);

  // Compute the flux at an interface
  void interface_flux_ (enzo_float *prim, enzo_float *cons, enzo_float *fluxes,
			enzo_float mag_pressure);
};

#endif /* ENZO_ENZO_RIEMANN_HPP */
