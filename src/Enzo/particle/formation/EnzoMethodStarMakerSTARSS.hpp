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
  EnzoMethodStarMakerSTARSS(ParameterGroup p);

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

  /// return minimum AMR level where star formation can occur
  ///
  /// @note
  /// The override keyword tells the compiler to perform a sanity-check for us
  /// compilation fails if we aren't overwriting a virtual function
  int sf_minimum_level() const noexcept override
  { return min_level_; }

protected:

  /// minimum AMR level for star formation
  int min_level_;
  bool turn_off_probability_;

};

#endif /* EnzoMethodStarMakerSTARSS */
