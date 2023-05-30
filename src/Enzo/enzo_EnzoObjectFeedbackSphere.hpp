class EnzoObjectFeedbackSphere : public ObjectSphere {
  // extension of Cello's ObjectSphere class to include metal yields
  public:
    /// Constructor
    EnzoObjectFeedbackSphere(double center[3], double radius, 
      double metal_yield_SNe, double metal_yield_HNe, double metal_yield_PISNe) 
      throw()
      : ObjectSphere(center, radius),
        metal_yield_SNe_(metal_yield_SNe),
        metal_yield_HNe_(metal_yield_HNe),
        metal_yield_PISNe_(metal_yield_PISNe)
        { };

    EnzoObjectFeedbackSphere() throw()
      : ObjectSphere()
    { };

    /// Charm++ PUP::able declarations
    PUPable_decl(EnzoObjectFeedbackSphere);

    EnzoObjectFeedbackSphere (CkMigrateMessage *m)
      : ObjectSphere(m),
        metal_yield_SNe_(0),
        metal_yield_HNe_(0),
        metal_yield_PISNe_(0)
    { }

    void pup (PUP::er &p)
    {
      TRACEPUP;

      ObjectSphere::pup(p);
      p | metal_yield_SNe_;
      p | metal_yield_HNe_;
      p | metal_yield_PISNe_; 
    }

  public:

    double metal_mass_SNe()   {return metal_yield_SNe_;}
    double metal_mass_HNe()   {return metal_yield_HNe_;}
    double metal_mass_PISNe() {return metal_yield_PISNe_;}

  private:

    double metal_yield_SNe_, metal_yield_HNe_, metal_yield_PISNe_;
};


