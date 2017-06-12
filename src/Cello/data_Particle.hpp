// See LICENSE_CELLO file for license and copyright information

/// @file     data_Particle.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-10-13
/// @brief    [\ref Data] Declaration of the Particle class
///
/// The Particle class is used to unify the interface of the global
/// ParticleDescr object and a Block's ParticleData object.

#ifndef DATA_PARTICLE_HPP
#define DATA_PARTICLE_HPP

class Particle {

  /// @class    Particle
  /// @ingroup  Data
  /// @brief    [\ref Data] 

public: // interface

  /// Constructor
  Particle() 
    : particle_descr_ (NULL),
      particle_data_ (NULL)
  { }

  /// Constructor
  Particle(ParticleDescr * particle_descr,
	   ParticleData  * particle_data) throw()
    : particle_descr_ (particle_descr),
      particle_data_ (particle_data)
  {
    if (particle_data_) particle_data_->allocate(particle_descr);
  }

  /// Copy constructor
  Particle(const Particle & particle) throw()
  {
    particle_descr_ = particle.particle_descr_;
    particle_data_ = particle.particle_data_; 
  }

  /// Assignment operator
  Particle & operator= (const Particle & particle) throw()
  { 
    particle_descr_ = particle.particle_descr_;
    particle_data_ = particle.particle_data_;
    return *this;
  }

  /// Comparison operator
  bool operator== (const Particle & particle) throw ()
  {
    return (particle_descr_ == particle.particle_descr_) &&
      (*particle_data_ == *particle.particle_data_);
  }

  /// Destructor
  ~Particle() throw()
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p)
  {
    WARNING ("Particle::pup()",
	     "Skipping since Particle is intended as transient object");
  };
  
  /// Return the particle descriptor for this particle
  ParticleDescr * particle_descr() { return particle_descr_; }

  /// Return the particle data for this particle
  ParticleData * particle_data() { return particle_data_; }

  //==================================================
  // ParticleDescr
  //==================================================

  /// Create a new type and return its id

  int new_type(std::string type)
  { 
    int it = particle_descr_->new_type(type); 
    particle_data_->allocate(particle_descr_);
    return it;
  }

  /// Return the number of types of particles

  int num_types() const
  { return particle_descr_->num_types(); }
  
  /// Return the index for the given particle type

  int type_index (std::string type) const
  { return particle_descr_->type_index(type); }

  /// Return the name of the given particle type given its index

  std::string type_name (int index) const
  { return particle_descr_->type_name(index); }

  //--------------------------------------------------
  // ATTRIBUTES
  //--------------------------------------------------

  /// Create a new attribute for the given type and return its id

  int new_attribute(int it, std::string attribute, int attribute_type)
  { return particle_descr_->new_attribute (it,attribute,attribute_type); }

  /// Return the number of attributes of the given type.

  int num_attributes(int it) const
  { return particle_descr_->num_attributes(it); }

  /// Return the index for the given attribute

  int attribute_index (int it, std::string attribute) const
  { return particle_descr_->attribute_index(it,attribute); }

  /// Return the name of the given attribute

  std::string attribute_name (int it, int ia) const
  { return particle_descr_->attribute_name(it,ia); }

  /// Define which attributes represent position coordinates (-1 if
  /// not defined)

  void set_position (int it, int ix, int iy=-1, int iz=-1)
  { particle_descr_->set_position (it,ix,iy,iz); }

  /// Define which attributes represent velocity coordinates (-1 if
  /// not defined)

  void set_velocity (int it, int ix, int iy=-1, int iz=-1)
  { particle_descr_->set_velocity (it,ix,iy,iz); }

  /// Byte offsets of attributes into block array.  Not including
  /// initial offset for 16-byte alignment.

  int attribute_offset(int it, int ia) const
  { return particle_descr_->attribute_offset(it,ia); }

  //--------------------------------------------------
  // CONSTANTS
  //--------------------------------------------------

  /// Create a new constant for the given type and return its id

  int new_constant(int it, std::string constant, int constant_type)
  { return particle_descr_->new_constant (it,constant,constant_type); }

  /// Return the number of constants of the given type.

  int num_constants(int it) const
  { return particle_descr_->num_constants(it); }

  /// Return the index for the given constant

  int constant_index (int it, std::string constant) const
  { return particle_descr_->constant_index(it,constant); }

  /// Return the name of the given constant

  std::string constant_name (int it, int ia) const
  { return particle_descr_->constant_name(it,ia); }

  /// Return the number of bytes allocated for the given constant.

  int constant_bytes (int it,int ic) const
  { return particle_descr_->constant_bytes(it,ic); }

  /// Return the constant array for the given particle type

  char * constant_array (int it)
  { return particle_descr_->constant_array(it); }

  /// Return a pointer to the given constant for the given type

  char * constant_value (int it, int ic)
  { return particle_descr_->constant_value(it, ic); }

  /// Byte offsets of constants into block array.

  int constant_offset(int it, int ic) const
  { return particle_descr_->constant_offset(it,ic); }


  //--------------------------------------------------
  // BYTES
  //--------------------------------------------------

  /// Return the data type of the given attribute.
  int attribute_type (int it,int ia) const
  { return particle_descr_->attribute_type(it,ia); }

  /// Return the number of bytes allocated for the given attribute.
  int attribute_bytes (int it,int ia) const
  { return particle_descr_->attribute_bytes(it,ia); }

  /// Return the number of bytes use to represent a particle.
  int particle_bytes (int it) const
  { return particle_descr_->particle_bytes(it); }

  /// Return the attribute corresponding to the given position
  /// coordinate, -1 if none
  int attribute_position (int it, int axis)
  { return particle_descr_->attribute_position_[it][axis]; }

  /// Return the attribute corresponding to the given velocity
  /// coordinate, -1 if none
  int attribute_velocity (int it, int axis)
  { return particle_descr_->attribute_velocity_[it][axis]; }

  /// Return the stride of the given attribute if interleaved, otherwise 1.
  /// Computed as attribute\_bytes(it) / attribute\_bytes(it,ia).
  /// Must be evenly divisible.

  int stride(int it, int ia) const
  { return particle_descr_->stride(it,ia); }

  //--------------------------------------------------
  // INTERLEAVING
  //--------------------------------------------------

  /// Set whether attributes are interleaved for the given type.

  void set_interleaved (int it, bool interleaved)
  { particle_descr_->set_interleaved(it,interleaved); }

  /// Return whether attributes are interleaved or not

  bool interleaved (int it) const
  { return particle_descr_->interleaved(it); }

  /// Return the current batch size.

  int batch_size() const
  { return particle_descr_->batch_size(); }

  /// Return the batch and particle index given a global particle
  /// index i.  This is useful e.g. for iterating over range of
  /// particles, e.g. initializing new particles after insert().
  /// Basically just div / mod.  ASSUMES COMPRESSED.

  void index (int i, int * ib, int * ip) const
  { particle_descr_->index(i,ib,ip); }

  /// Return the Grouping object for the particle types

  Grouping * groups ()
  { return particle_descr_->groups(); }

  //==================================================
  // ParticleData
  //==================================================

  /// Return the attribute array for the given particle type and batch

  char * attribute_array (int it,int ia,int ib)
  { return particle_data_->attribute_array 
      (particle_descr_, it,ia,ib); }

  /// Return the number of batches of particles for the given type.

  int num_batches (int it) const
  { return particle_data_->num_batches(it); }

  /// Return the number of particles in the given batch, of the given
  /// type, or total on the block.

  int num_particles (int it, int ib) const
  { return particle_data_->num_particles(particle_descr_,it,ib); }
  int num_particles (int it) const
  { return particle_data_->num_particles(particle_descr_,it); }
  int num_particles () const
  { return particle_data_->num_particles(particle_descr_); }

  /// Create the given number of particles of the given type.  Always
  /// creates them at the end instead of filling up any unused
  /// particle spaces in earlier batches, to ease initialization via
  /// index()

  int insert_particles (int it, int np)
  { return  particle_data_->insert_particles (particle_descr_, it, np); }

  /// Delete the given particles in the batch according to mask
  /// attribute.  Compresses the batch if particles deleted, so batch
  /// may have fewer than max number of particles.  Other batches
  /// remain unchanged.  Return the number of particles deleted.

  int delete_particles (int it, int ib, const bool * mask = 0)
  { return particle_data_->delete_particles (particle_descr_,it,ib,mask); }

  /// Scatter particles among an array of other Particle structures.
  /// Typically used for preparing to send particles that have gone
  /// out of the block to neighboring blocks.

  void scatter (int it, int ib,
		int np, const bool * mask, const int * index,
		int n, ParticleData ** particle_array)
  { particle_data_->scatter
      (particle_descr_,it,ib,np,mask,index,n,particle_array);  }
  
  /// Gather particles from an array of other Particle structures.
  /// Typically used after receiving particles from neighboring blocks
  /// that have entered this block.  Return the total number of particles
  /// inserted

  int gather (int it, int n, ParticleData **particle_array)
  { return particle_data_->gather(particle_descr_,it,n,particle_array); }

  /// Compress particles in batches so that all batches except
  /// possibly the last have batch_size() particles.  May be performed
  /// periodically to recover unused memory from multiple insert/deletes

  void compress ()
  { particle_data_->compress(particle_descr_); }
  void compress (int it)
  { particle_data_->compress(particle_descr_,it); }

  /// Return the storage "efficiency" for particles of the given type
  /// and in the given batch, or average if batch or type not specified.
  /// 1.0 means no wasted storage, 0.5 means twice as much storage
  /// is being used.  Defined as 1 / overhead().

  float efficiency ()
  { return particle_data_->efficiency(particle_descr_); }
  float efficiency (int it)
  { return particle_data_->efficiency(particle_descr_,it); }
  float efficiency (int it, int ib)
  { return particle_data_->efficiency(particle_descr_,it,ib); }

  /// Return the storage "overhead" for particles of the given type
  /// and in the given batch, or average if batch or type not specified.
  /// 1.0 means no wasted storage, 2.0 means twice as much storage
  /// is being used.  Defined as 1 / overhead().

  float overhead ()
  { return particle_data_->overhead(particle_descr_); }
  float overhead (int it)
  { return particle_data_->overhead(particle_descr_,it); }
  float overhead (int it, int ib)
  { return particle_data_->overhead(particle_descr_,it,ib); }

  /// Fill a vector of position coordinates for the given type and batch
  bool position (int it, int ib,
		 double * x, double * y = 0, double * z = 0)
  { return particle_data_->position(particle_descr_,it,ib,x,y,z); }

  void position_update (int it, int ib, 
			long double dx, long double dy, long double dz)
  { particle_data_->position_update (particle_descr_,it,ib,dx,dy,dz);  }


  /// Fill a vector of velocity coordinates for the given type and batch
  bool velocity (int it, int ib,
		 double * vx, double * vy = 0, double * vz = 0)
  { return particle_data_->velocity(particle_descr_,it,ib,vx,vy,vz); }

  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size () const
  { return particle_data_->data_size (particle_descr_); }

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (char * buffer) const
  { return particle_data_->save_data (particle_descr_,buffer); }

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (char * buffer)
  { return particle_data_->load_data (particle_descr_,buffer); }

  //--------------------------------------------------

  void debug () const
  { particle_data_->debug (particle_descr_); }

  void write_ifrite (int it, std::string file_name,
		     double xm, double ym, double zm,
		     double xp, double yp, double zp) 
  { particle_data_->write_ifrite 
      (particle_descr_,it,file_name,xm,ym,zm,xp,yp,zp);  }

private: // functions

  /// Return an id (not "index"); for a particle that is guaranteed to
  /// be unique across all processors.  May involve communication.

  /// long long assign_id_ ()

private: // attributes

  /// Particle descriptor for global particle data
  ParticleDescr * particle_descr_;

  /// Particle data for the specific Block
  ParticleData * particle_data_;

  // NOTE: change pup() function whenever attributes change

};

#endif /* DATA_PARTICLE_HPP */
