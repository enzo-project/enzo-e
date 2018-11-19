
#ifndef ENZO_ENZO_RIEMANN_HPP
#define ENZO_ENZO_RIEMANN_HPP
class EnzoRiemann : public PUP::able
{

  /// @class    EnzoRiemann
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate approximate Riemann Solvers

public: // interface

  /// Create a new EnzoRiemann object
  EnzoRiemann() throw()
  {}

  /// CHARM++ PUP::able declaration
  //PUPable_abstract(EnzoRiemann);

  /// CHARM++ migration constructor for PUP::able
  //EnzoEquationOfState (CkMigrateMessage *m)
  //  : PUP::able(m)
  //{  }

  /// Virtual destructor
  virtual ~EnzoRiemann()
  {  }

  /// CHARM++ Pack / Unpack function
  //void pup (PUP::er &p);

  /// Solve the Riemann Problem
  /// Probably want to make dim, an enum (x,y or z) Current plan: [0, 1, 2]
  /// dim tells the solver which dimension to compute fluxes and indicates the
  /// dimension along which the reconstructed primitive values are face-centered
  virtual void solve (Block *block, std::vector<int> &priml_ids,
		      std::vector<int> &primr_ids,
		      std::vector<int> &flux_ids, int dim,
		      EnzoEquationOfState *eos)=0;
};

#endif /* ENZO_ENZO_RIEMANN_HPP */
