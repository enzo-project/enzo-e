
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
  //EnzoRiemann (CkMigrateMessage *m)
  //  : PUP::able(m)
  //{  }

  /// Virtual destructor
  virtual ~EnzoRiemann()
  {  }

  /// CHARM++ Pack / Unpack function
  //void pup (PUP::er &p);

  /// Solve the Riemann Problem - dim (0, 1, or 2) tells the solver which
  /// dimension to compute fluxes along and indicates the dimension along which
  /// the reconstructed primitive values are face-centered
  virtual void solve (Block *block, Grouping &priml_group,
		      Grouping &primr_group, Grouping &flux_group, int dim,
		      EnzoEquationOfState *eos)=0;
};

#endif /* ENZO_ENZO_RIEMANN_HPP */
