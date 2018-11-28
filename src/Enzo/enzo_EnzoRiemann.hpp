
#ifndef ENZO_ENZO_RIEMANN_HPP
#define ENZO_ENZO_RIEMANN_HPP

// All Riemann solvers are expected to calculate fluxes at all cell faces that
// that are not on the outermost edges. If there are N cells along the
// direction of flux calculations, the only 2 faces where flux is not computed
// is at i=-1/2 and i=N+1/2 (the index of the first cell-center is i=0)

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
