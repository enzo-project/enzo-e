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
  Particle(ParticleDescr * particle_descr,
	   ParticleData  * particle_data) throw()
    : particle_descr_ (particle_descr),
      particle_data_ (particle_data)
  {}

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

  /// Destructor
  ~Particle() throw()
  {};

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
  { return particle_descr_->new_type(type); }

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

  //--------------------------------------------------
  // BYTES
  //--------------------------------------------------

  /// Return the number of bytes per particle allocated for all attributes

  int attribute_bytes (int it) const
  { return particle_descr_->attribute_bytes(it); }

  /// Return the number of bytes allocated for the given attribute.

  int attribute_bytes(int it, int ia) const
  { return particle_descr_->attribute_bytes(it,ia); }

  //--------------------------------------------------
  // INTERLEAVING
  //--------------------------------------------------

  /// Set whether attributes are interleaved for the given type.

  void set_interleaved (int it, bool interleaved)
  { particle_descr_->set_interleaved(it,interleaved); }

  /// Return whether attributes are interleaved or not

  bool interleaved (int it) const
  { return particle_descr_->interleaved(it); }

  /// Return the stride of the given attribute if interleaved, otherwise 1.
  /// Computed as attribute\_bytes(it) / attribute\_bytes(it,ia).
  /// Must be evenly divisible.

  int stride(int it, int ia) const
  { return particle_descr_->stride(it,ia); }

  /// Set the size of batches.  Must be set at most once.  May be
  /// defined when ParticleDescr is created.

  // int set_batch_size (mb)

  /// Return the batch and particle index given a global particle
  /// index i.  (useful for iterating over range of particles,
  /// e.g. initializing new particles after insert(); Basically just
  /// div / mod.

  // int index (i, it, &ib, &ip)

  /// Return the Grouping object for the particle types

  // Grouping * groups ();

  //==================================================
  // ParticleData
  //==================================================

  /// Return the attribute array for the given particle type and batch

  /// char * attribute_array (it,ia,ib)

  /// Return the number of batches of particles for the given type.

  /// int num_batches (it)

  /// Return the number of particles in the given batch, of the given
  /// type, or total on the block.

  /// int num_particles (it,ib)
  /// int num_particles (it)
  /// int num_particles ()

  /// Create the given number of particles of the given type.  Always
  /// creates them at the end instead of filling up any unused
  /// particle spaces in earlier batches, to ease initialization via
  /// index()

  /// void insert_particles (it,np)

  /// Delete the given particles in the batch according to mask
  /// attribute.  Compresses the batch if particles deleted, so batch
  /// may have fewer than max number of particles.  Other batches
  /// remain unchanged.

  /// void delete_particles (it,ib,im)

  /// Same as delete, but inserts particles into a second Particle
  /// object.  void split_particles (it,ib,im,particle) Compress
  /// particles in batches so that ib'th batch size equals batch_size.
  /// May be performed periodically to recover space lost in multiple
  /// insert/deletes

  /// void compress (it)

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
