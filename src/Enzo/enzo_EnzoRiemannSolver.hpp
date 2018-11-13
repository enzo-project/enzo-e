
#ifndef ENZO_ENZO_RIEMANNSOLVER_HPP
#define ENZO_ENZO_RIEMANNSOLVER_HPP
class EnzoRiemannSolver : public PUP::able
{

  /// @class    EnzoRiemannSolver
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Encapsulate Approximate Riemann Solvers

public: // interface

  /// Create a new EnzoRiemannSolver object
  EnzoRiemannSolver() throw()
  {}

  /// CHARM++ PUP::able declaration
  //PUPable_abstract(EnzoRiemannSolver);

  /// CHARM++ migration constructor for PUP::able
  //EnzoEquationOfState (CkMigrateMessage *m)
  //  : PUP::able(m)
  //{  }

  /// Virtual destructor
  virtual ~EnzoRiemannSolver()
  {  }

  /// CHARM++ Pack / Unpack function
  //void pup (PUP::er &p);

  /// Solve the Riemann Problem
  /// Probably want to make dim, an enum (x,y or z) Current plan: [0, 1, 2]
  /// dim tells the solver which dimension to compute fluxes and indicates the
  /// dimension along which the reconstructed primitive values are face-centered
  virtual void solve (Block *block, const std::vector<int> &priml_ids,
		      const std::vector<int> &primr_ids,
		      std::vector<int> &flux_ids, int dim,
		      EnzoEquationOfState *eos)=0;
};

#endif /* ENZO_ENZO_RIEMANNSOLVER_HPP */
