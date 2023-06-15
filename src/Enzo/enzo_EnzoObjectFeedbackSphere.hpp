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
      : ObjectSphere(),
        metal_yield_SNe_(0),
        metal_yield_HNe_(0),
        metal_yield_PISNe_(0)
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

    return pc;
  }

  //---------------------------------------------------------------------
 
  public:

    double metal_mass_SNe()   {return metal_yield_SNe_;}
    double metal_mass_HNe()   {return metal_yield_HNe_;}
    double metal_mass_PISNe() {return metal_yield_PISNe_;}

  private:

    double metal_yield_SNe_, metal_yield_HNe_, metal_yield_PISNe_;
};


