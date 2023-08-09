class EnzoObjectFeedbackSphere : public ObjectSphere {
  // extension of Cello's ObjectSphere class to include metal yields
  public:
    /// Constructor
    EnzoObjectFeedbackSphere(double center[3], double radius, 
      double metal_yield_SNe, double metal_yield_HNe, double metal_yield_PISNe,
      int num_SNe, int num_HNe, int num_PISNe, int num_BH, double stellar_mass) 
      throw()
      : ObjectSphere(center, radius),
        metal_yield_SNe_(metal_yield_SNe),
        metal_yield_HNe_(metal_yield_HNe),
        metal_yield_PISNe_(metal_yield_PISNe),
        num_SNe_(num_SNe),
        num_HNe_(num_HNe),
        num_PISNe_(num_PISNe),
        num_BH_(num_BH),
        stellar_mass_(stellar_mass)
        { };

    EnzoObjectFeedbackSphere() throw()
      : ObjectSphere(),
        metal_yield_SNe_(0),
        metal_yield_HNe_(0),
        metal_yield_PISNe_(0),
        num_SNe_(0),
        num_HNe_(0),
        num_PISNe_(0),
        num_BH_(0),
        stellar_mass_(0)
    { };

    /// Charm++ PUP::able declarations
    PUPable_decl(EnzoObjectFeedbackSphere);

    EnzoObjectFeedbackSphere (CkMigrateMessage *m)
      : ObjectSphere(m),
        metal_yield_SNe_(0),
        metal_yield_HNe_(0),
        metal_yield_PISNe_(0),
        num_SNe_(0),
        num_HNe_(0),
        num_PISNe_(0),
        num_BH_(0),
        stellar_mass_(0)
    { }

    void pup (PUP::er &p)
    {
      TRACEPUP;

      ObjectSphere::pup(p);
      p | metal_yield_SNe_;
      p | metal_yield_HNe_;
      p | metal_yield_PISNe_;
      p | num_SNe_;
      p | num_HNe_;
      p | num_PISNe_;
      p | num_BH_;
      p | stellar_mass_; 
    }

  ///--------------------
  /// PACKING / UNPACKING
  ///--------------------
  
  /// Return the number of bytes required to serialize the data object
  int data_size ()
  {
    //--------------------------------------------------
    //  1. determine buffer size (must be consistent with #3)
    //--------------------------------------------------

    int size = 0;

    double center_[3] = {this->pos(0),this->pos(1),this->pos(2)};
    double radius_ = this->r();
    SIZE_ARRAY_TYPE(size,double,center_,3);
    SIZE_SCALAR_TYPE(size,double,radius_);
    SIZE_SCALAR_TYPE(size,double,metal_yield_SNe_);
    SIZE_SCALAR_TYPE(size,double,metal_yield_HNe_);
    SIZE_SCALAR_TYPE(size,double,metal_yield_PISNe_);
    SIZE_SCALAR_TYPE(size,int,num_SNe_);
    SIZE_SCALAR_TYPE(size,int,num_HNe_);
    SIZE_SCALAR_TYPE(size,int,num_PISNe_);
    SIZE_SCALAR_TYPE(size,int,num_BH_);
    SIZE_SCALAR_TYPE(size,double,stellar_mass_);

    return size;
  }

  //----------------------------------------------------------------------

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer)
  {
    char * pc = buffer;
    double center_[3] = {this->pos(0),this->pos(1),this->pos(2)};
    double radius_ = this->r();
    SAVE_ARRAY_TYPE(pc,double,center_,3);
    SAVE_SCALAR_TYPE(pc,double,radius_);
    SAVE_SCALAR_TYPE(pc,double,metal_yield_SNe_);
    SAVE_SCALAR_TYPE(pc,double,metal_yield_HNe_);
    SAVE_SCALAR_TYPE(pc,double,metal_yield_PISNe_);
    SAVE_SCALAR_TYPE(pc,int,num_SNe_);
    SAVE_SCALAR_TYPE(pc,int,num_HNe_);
    SAVE_SCALAR_TYPE(pc,int,num_PISNe_);
    SAVE_SCALAR_TYPE(pc,int,num_BH_);
    SAVE_SCALAR_TYPE(pc,double,stellar_mass_);

    ASSERT2 ("EnzoObjectFeedbackSphere::save_data()",
             "Expecting buffer size %d actual size %d",
             data_size(),(pc-buffer),
             (data_size() == (pc-buffer)));

    return pc;
  }

  //----------------------------------------------------------------------

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer)
  {
    char * pc = buffer;
    double center_[3] = {this->pos(0),this->pos(1),this->pos(2)};
    double radius_ = this->r();
    LOAD_ARRAY_TYPE(pc,double,center_,3);
    LOAD_SCALAR_TYPE(pc,double,radius_);
    LOAD_SCALAR_TYPE(pc,double,metal_yield_SNe_);
    LOAD_SCALAR_TYPE(pc,double,metal_yield_HNe_);
    LOAD_SCALAR_TYPE(pc,double,metal_yield_PISNe_);
    LOAD_SCALAR_TYPE(pc,int,num_SNe_);
    LOAD_SCALAR_TYPE(pc,int,num_HNe_);
    LOAD_SCALAR_TYPE(pc,int,num_PISNe_);
    LOAD_SCALAR_TYPE(pc,int,num_BH_);
    LOAD_SCALAR_TYPE(pc,double,stellar_mass_);

    return pc;
  }

  //---------------------------------------------------------------------
 
  public:

    double metal_mass_SNe()   {return metal_yield_SNe_;}
    double metal_mass_HNe()   {return metal_yield_HNe_;}
    double metal_mass_PISNe() {return metal_yield_PISNe_;}
    double stellar_mass()     {return stellar_mass_;}
    int num_SNe()   {return num_SNe_;}
    int num_HNe()   {return num_HNe_;}
    int num_PISNe() {return num_PISNe_;}
    int num_BH()    {return num_BH_;}
  private:

    double metal_yield_SNe_, metal_yield_HNe_, metal_yield_PISNe_;
    double stellar_mass_;
    int num_SNe_, num_HNe_, num_PISNe_, num_BH_;
};


