// See LICENSE_CELLO file for license and copyright information

/// @file     data_Data.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2013-03-10
/// @brief    [\ref Data] Declaration of the Data class
///

#ifndef DATA_DATA_HPP
#define DATA_DATA_HPP

class Data {

  friend class Block;
  friend class IoBlock;

  /// @class    Data
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Constructor
  Data(int nx, int ny, int nz,
       int num_field_data,
       double xm, double xp,
       double ym, double yp,
       double zm, double zp,
       // specifying field_descr and particle_descr included for test code only
       FieldDescr * field_descr = nullptr, 
       ParticleDescr * particle_descr = nullptr) throw();

  /// Destructor
  ~Data() throw();

  /// Copy constructor
  Data(const Data & data) throw();

  /// Assignment operator
  Data & operator= (const Data & data) throw();

  /// Empty constructor
  Data()
    : num_field_data_(0),
      field_data_(),
      particle_data_(nullptr),
      flux_data_(nullptr),
      scalar_data_long_double_(),
      scalar_data_double_(),
      scalar_data_int_(),
      scalar_data_long_long_(),
      scalar_data_sync_(),
      scalar_data_void_(),
      scalar_data_index_()
  {
    lower_[0] = 0.0;
    lower_[1] = 0.0;
    lower_[2] = 0.0;
    upper_[0] = 0.0;
    upper_[1] = 0.0;
    upper_[2] = 0.0;
  }

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    bool up = p.isUnpacking();
    p | num_field_data_;
    // allocate field_data_[] vector first if unpacking
    if (up) field_data_.resize(num_field_data_);
    for (int i=0; i<num_field_data_; i++) {
      if (up) field_data_[i] = new FieldData;
      p | *field_data_[i];
    }
    if (up) {
      particle_data_ = new ParticleData;
    }
    p | *particle_data_;
    if (up) {
      flux_data_ = new FluxData;
    }
    p | *flux_data_;
    p | scalar_data_long_double_;
    p | scalar_data_double_;
    p | scalar_data_int_;
    p | scalar_data_long_long_;
    p | scalar_data_sync_;
    //    p | scalar_data_void_;
    p | scalar_data_index_;
    static bool warn[CONFIG_NODE_SIZE] = {false};
    const int in = cello::index_static();
    if (! warn[in]) {
      WARNING("Data::pup()","Skipping scalar_data_void_");
      warn[in]=true;
    }
    PUParray(p,lower_,3);
    PUParray(p,upper_,3);
    // NOTE: change this function whenever attributes change
  }

  /// Return domain lower extent
  inline void lower(double * x, 
		    double * y = 0,
		    double * z = 0) const throw ()
  {
    if (x) *x = lower_[0];
    if (y) *y = lower_[1];
    if (z) *z = lower_[2];
  }

  //----------------------------------------------------------------------

  /// Return domain upper extent
  inline void upper(double * x,
		    double * y = 0,
		    double * z = 0) const throw ()
  {
    if (x) *x = upper_[0];
    if (y) *y = upper_[1];
    if (z) *z = upper_[2];
  }

  //----------------------------------------------------------------------

  void allocate () throw();

  //----------------------------------------------------------------------
  // fields
  //----------------------------------------------------------------------

  /// Return the number of FieldData
  int num_field_data() const throw()
  { return num_field_data_; }

  /// Return the ith Field data
  const FieldData * field_data (size_t i=0) const throw()
  { return (i < field_data_.size()) ? field_data_[i] : NULL; }

  /// Return the ith Field data
  FieldData * field_data (size_t i=0) throw()
  { return (i < field_data_.size()) ? field_data_[i] : NULL; }

  /// Return the ith Field
  Field field (size_t i=0) throw()
  { return Field(cello::field_descr(),field_data(i)); }

  /// Return the x,y,z,t coordinates of field cell centers
  void field_cells (double * x, double * y, double * z,
		    int gx = 0, int gy = 0, int gz = 0) const
  { return field_cell_faces(x,y,z,gx,gy,gz,0,0,0); }

  /// Return the x,y,z,t coordinates of field cell faces (if there are n cell
  /// centers, then there are n+1 cell faces)
  void field_cell_faces (double * x, double * y, double * z,
			 int gx = 0, int gy = 0, int gz = 0,
			 int cx = 1, int cy = 1, int cz = 1) const;

  /// Return the cell widths of Fields
  void field_cell_width (double * hx, 
			 double * hy = 0,
			 double * hz = 0) const;

  //----------------------------------------------------------------------
  // particles
  //----------------------------------------------------------------------

  /// Return the constant Particle data
  const ParticleData * particle_data () const throw()
  { return particle_data_; }

  /// Return the Particle data
  ParticleData * particle_data () throw()
  { return particle_data_; }

  // Return the ith Particle descriptor
  ParticleDescr * particle_descr () throw();
  const ParticleDescr * particle_descr () const throw();

  /// Return the Particle object
  Particle particle () throw()
  { return Particle(cello::particle_descr(),
		    particle_data_); }

  //----------------------------------------------------------------------
  // fluxes
  //----------------------------------------------------------------------

  /// Return the constant Flux data
  const FluxData * flux_data () const throw()
  { return flux_data_; }

  /// Return the Flux data
  FluxData * flux_data () throw()
  { return flux_data_; }

  //----------------------------------------------------------------------
  // scalars
  //----------------------------------------------------------------------

  /// Return the scalar_data objects
  
  ScalarData<double> * scalar_data_double ()
  { return &scalar_data_double_; }
  ScalarData<long double> * scalar_data_long_double ()
  { return &scalar_data_long_double_; }
  ScalarData<int> * scalar_data_int ()
  { return &scalar_data_int_; }
  ScalarData<long long> * scalar_data_long_long ()
  { return &scalar_data_long_long_; }
  ScalarData<Sync> * scalar_data_sync ()
  { return &scalar_data_sync_; }
  ScalarData<void *> * scalar_data_void ()
  { return &scalar_data_void_; }
  ScalarData<Index> * scalar_data_index ()
  { return &scalar_data_index_; }

  /// Return the Scalar objects
  Scalar<long double> scalar_long_double()
  { return Scalar<long double>
      (cello::scalar_descr_long_double(),
       &scalar_data_long_double_); }
  Scalar<double> scalar_double()
  { return Scalar<double>
      (cello::scalar_descr_double(),
       &scalar_data_double_); }
  Scalar<int> scalar_int()
  { return Scalar<int>
      (cello::scalar_descr_int(),
       &scalar_data_int_); }
  Scalar<long long> scalar_long_long()
  { return Scalar<long long>
      (cello::scalar_descr_long_long(),
       &scalar_data_long_long_); }
  Scalar<Sync> scalar_sync()
  { return Scalar<Sync>
      (cello::scalar_descr_sync(),
       &scalar_data_sync_); }
  Scalar<void *> scalar_void()
  { return Scalar<void *>
      (cello::scalar_descr_void (),
       &scalar_data_void_); }
  Scalar<Index> scalar_index()
  { return Scalar<Index>
      (cello::scalar_descr_index(),
       &scalar_data_index_); }


private: // functions

  void copy_(const Data & data) throw();

private: // attributes

  /// Number of field_data's (required by CHARM++ PUP::er)
  int num_field_data_;

  /// Array of field data
  std::vector<FieldData *> field_data_;

  /// Particle data
  ParticleData * particle_data_;

  /// Flux data
  FluxData * flux_data_;

  /// Scalar data
  ScalarData<long double> scalar_data_long_double_;
  ScalarData<double>      scalar_data_double_;
  ScalarData<int>         scalar_data_int_;
  ScalarData<long long>   scalar_data_long_long_;
  ScalarData<Sync>        scalar_data_sync_;
  ScalarData<void *>      scalar_data_void_;
  ScalarData<Index>       scalar_data_index_;

  /// Lower extent of the box associated with the block [computable]
  double lower_[3];

  /// Upper extent of the box associated with the block [computable]
  double upper_[3];

  // NOTE: change pup() function whenever attributes change

};

#endif /* DATA_DATA_HPP */

