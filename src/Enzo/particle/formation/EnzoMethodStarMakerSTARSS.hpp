///

///
///
///
///

#ifndef ENZO_ENZO_METHOD_STARMAKER_STARSS
#define ENZO_ENZO_METHOD_STARMAKER_STARSS

///
///
///
///

class EnzoMethodStarMakerSTARSS : public EnzoMethodStarMaker {

public:
  // Create new EnzoStarMakerSTARSS object
  EnzoMethodStarMakerSTARSS();

  // Charm++ PUP::able declarations
  PUPable_decl(EnzoMethodStarMakerSTARSS);

  // Charm++ PUP::able declarations
  EnzoMethodStarMakerSTARSS (CkMigrateMessage *m)
   : EnzoMethodStarMaker (m)
   {  }

  /// Charm++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Apply method
  virtual void compute ( Block * block) throw();

  virtual std::string particle_type () throw()
  { return "star";}

  /// Name
  virtual std::string name () throw()
   { return "star_maker";}

  virtual ~EnzoMethodStarMakerSTARSS() throw() {};

protected:

};

#endif /* EnzoMethodStarMakerSTARSS */
