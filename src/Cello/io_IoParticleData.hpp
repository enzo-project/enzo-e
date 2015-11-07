// See LICENSE_CELLO file for license and copyright information

/// @file     io_IoParticleData.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2011-10-06
/// @brief    [\ref Io] Declaration of the IoParticleData class

#ifndef IO_IO_PARTICLE_DATA_HPP
#define IO_IO_PARTICLE_DATA_HPP

class ParticleData;
class ParticleDescr;

class IoParticleData : public Io {

  /// @class    IoParticleData
  /// @ingroup  Io
  /// @brief    [\ref Io] Class for linking between ParticleData and Output classes

public: // interface

  /// Constructor
  IoParticleData(const ParticleDescr * particle_descr) throw();

  /// Destructor
  virtual ~IoParticleData() throw()
  {}

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

  /// Set ParticleIndex
  void set_particle_index (int particle_index) throw()
  { particle_index_ = particle_index;};

  /// Set ParticleData
  void set_particle_data (ParticleData * particle_data) throw()
  { particle_data_ = particle_data;};

  /// Set ParticleDescr
  void set_particle_descr (ParticleDescr * particle_descr) throw()
  { particle_descr_ = particle_descr; };


#include "_io_Io_common.hpp"

  
protected: // functions

  /// ParticleDescr for the ParticleData
  ParticleDescr * particle_descr_;

  /// Current ParticleData
  ParticleData * particle_data_;

  /// Index of the particle in the ParticleData
  int particle_index_;


private: // attributes


};

#endif /* IO_IO_PARTICLE_DATA_HPP */

