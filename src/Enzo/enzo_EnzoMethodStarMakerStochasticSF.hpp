///

///
///
///
///

#ifndef ENZO_ENZO_METHOD_STARMAKER_STOCHASTICSF
#define ENZO_ENZO_METHOD_STARMAKER_STOCHASTICSF

///
///
///
///

class EnzoMethodStarMakerStochasticSF : public EnzoMethodStarMaker {

public:
  // Create new EnzoStarMakerStochasticSF object
  EnzoMethodStarMakerStochasticSF();

  // Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodStarMakerStochasticSF);

  // Charm++ PUP::able declarations
  EnzoMethodStarMakerStochasticSF (CkMigrateMessage *m)
   : EnzoMethodStarMaker (m)
   {  }

  /// Charm++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply method
  virtual void compute ( Block * block) throw();

  //virtual double timestep (Block * block) const throw();

  virtual std::string particle_type () throw()
  { return "star";}

  /// Name
  virtual std::string name () throw()
   { return "star_maker";}

  virtual ~EnzoMethodStarMakerStochasticSF() throw() {};

protected:

};

#endif /* EnzoMethodStarMakerStochasticSF */
