// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-10-12
/// @brief    [\ref Data] Declaration of the ParticleData class

#ifndef DATA_PARTICLE_DATA_HPP
#define DATA_PARTICLE_DATA_HPP

class ParticleData {

  /// @class    ParticleData
  /// @ingroup  Data
  /// @brief    [\ref Data] 

  friend class Particle;

public: // interface

  static int64_t counter[CONFIG_NODE_SIZE];

  /// Constructor
  ParticleData();

  /// Comparison operator
  bool operator== (const ParticleData & particle_data) throw ();

  /// Destructor
  ~ParticleData();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Return the attribute array for the given particle type and batch
  char * attribute_array (ParticleDescr *, int it, int ia, int ib);

  /// Return the number of batches of particles for the given type.

  int num_batches (int it) const;

  /// Return the number of particles in the given batch, of the given
  /// type, or total on the block.

  int num_particles (ParticleDescr *) const;
  int num_particles (ParticleDescr *, int it) const;
  int num_particles (ParticleDescr *, int it, int ib) const;

  /// Create the given number of particles of the given type.  Always
  /// creates them at the end instead of filling up any unused
  /// particle spaces in earlier batches, to ease initialization via
  /// index()

  int insert_particles (ParticleDescr *, int it, int np);

  /// Delete the given particles in the batch according to mask
  /// attribute.  Compresses the batch if particles deleted, so batch
  /// may have fewer than max number of particles.  Other batches
  /// remain unchanged.

  int delete_particles (ParticleDescr *, int it, int ib, const bool * m = NULL);

  /// Scatter particles among an array of other Particle structures.
  /// Typically used for preparing to send particles that have gone
  /// out of the block to neighboring blocks.  Particles are not deleted,
  /// so must be done so manually, e.g. using delete_particles.

  void scatter (ParticleDescr *, int it, int ib,
		int np, const bool * mask, const int * index,
		int n,  ParticleData * particle_array[]);
  
  /// Gather particles from an array of other Particle structures.
  /// Typically used after receiving particles from neighboring blocks
  /// that have entered this block.  Return the total number of particles
  /// inserted

  int gather (ParticleDescr *, int it, int n, ParticleData * particle_array[]);

  /// Compress particles in batches so that all batches except
  /// possibly the last have batch_size() particles.  May be performed
  /// periodically to recover unused memory from multiple insert/deletes

  void compress (ParticleDescr *);
  void compress (ParticleDescr *, int it);

  /// Return the storage "efficiency" for particles of the given type
  /// and in the given batch, or average if batch or type not specified.
  /// 1.0 means no wasted storage, 0.5 means twice as much storage
  /// is being used.  Defined as 1 / overhead().

  float efficiency (ParticleDescr *);
  float efficiency (ParticleDescr *, int it);
  float efficiency (ParticleDescr *, int it, int ib);

  /// Return the storage "overhead" for particles of the given type
  /// and in the given batch, or average if batch or type not specified.
  /// 1.0 means no wasted storage, 2.0 means twice as much storage
  /// is being used.  Defined as 1 / overhead().

  float overhead (ParticleDescr *);
  float overhead (ParticleDescr *, int it);
  float overhead (ParticleDescr *, int it, int ib);

  /// Add a new type to the attribute_array.  Should only be called
  /// by Particle.

  /// Allocate arrays based on ParticleDescr attributes
  void allocate(ParticleDescr * particle_descr)
  {
    const size_t nt = particle_descr->num_types();
    if (attribute_array_.size() < nt) {
      attribute_array_.resize(nt);
    }
    if (attribute_align_.size() < nt) {
      attribute_align_.resize(nt);
    }
    if (particle_count_.size() < nt) {
      particle_count_.resize(nt);
    }
  };

  /// Fill a vector of position coordinates for the given type and batch
  bool position (ParticleDescr * particle_descr,
		 int it, int ib,
		 double * x, double * y = 0, double * z = 0);

  /// Update positions given an increment
  void position_update 
  (ParticleDescr * particle_descr,int it, int ib, 
   long double dx, long double dy, long double dz);
			 

  /// Fill a vector of velocity coordinates for the given type and batch
  bool velocity (ParticleDescr * particle_descr,
		 int it, int ib,
		 double * vx, double * vy = 0, double * vz = 0);

  //--------------------------------------------------

  /// Return the number of bytes required to serialize the data object
  int data_size (ParticleDescr * particle_descr) const;

  /// Serialize the object into the provided empty memory buffer.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * save_data (ParticleDescr * particle_descr, char * buffer) const;

  /// Restore the object from the provided initialized memory buffer data.
  /// Returns the next open position in the buffer to simplify
  /// serializing multiple objects in one buffer.
  char * load_data (ParticleDescr * particle_descr, char * buffer);

  //--------------------------------------------------

  void debug (ParticleDescr * particle_descr);

private: /// functions

  /// Return an id (not "index"); for a particle that is guaranteed to
  /// be unique across all processors.  May involve communication.

  /// long long assign_id_ ()

  /// Allocate attribute_array_ block, aligned at 16 byte boundary
  /// with updated attribute_align_
  void resize_attribute_array_ (ParticleDescr *, int it, int ib, int np);

  void check_arrays_ (ParticleDescr * particle_descr,
		      std::string file, int line) const;

  /// Copy the given floating point attribute of given type (float,
  /// double, quad, etc.) to the given coordinate double position
  /// array.
  void copy_attribute_float_ 
  (ParticleDescr * particle_descr,
   int type, int it, int ib, int ia, double * coord);

  /// Increment the given floating point attribute of given float type
  /// (float, double, quad, etc.) by the given double long constant value.
  void update_attribute_float_ 
  (ParticleDescr * particle_descr,
   int type, int it, int ib, int ia, long double da);

  /// Positions may be defined relative to the Block using integer
  /// variables.  If so, copy the position coordinate to the given
  /// double position array, where coordinate range in the Block is
  /// [-1.0, 1.0).  Conversion to double may involve floating-point
  /// errors, especially if position is defined using large ints,
  /// e.g. int64_t
  void copy_position_int_ 
  (ParticleDescr * particle_descr,
   int type, int it, int ib, int ia, double * coord);

  /// Update version of copy_position_int_()
  void update_position_int_ 
  (ParticleDescr * particle_descr,
   int type, int it, int ib, int ia, int64_t da);

  void write_ifrite (ParticleDescr * particle_descr,
		     int it, std::string file_name,
		     double xm, double ym, double zm,
		     double xp, double yp, double zp);

private: /// attributes

  /// Array of blocks of particle attributes array_[it][ib][iap];
  std::vector< std::vector< std::vector<char> > > attribute_array_;

  /// Alignment adjustment to correct for 16-byte alignment of
  /// first attribute in each batch

  std::vector< std::vector< char > > attribute_align_;

  /// Number of particles in the batch particle_count_[it][ib];
  std::vector < std::vector < int > > particle_count_;

};

#endif /* DATA_PARTICLE_DATA_HPP */

