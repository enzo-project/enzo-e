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

public: // interface

  /// Constructor
  ParticleData();

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p);

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

private: /// functions

  /// Return an id (not "index"); for a particle that is guaranteed to
  /// be unique across all processors.  May involve communication.

  /// long long assign_id_ ()

private: /// attributes

};

#endif /* DATA_PARTICLE_DATA_HPP */

