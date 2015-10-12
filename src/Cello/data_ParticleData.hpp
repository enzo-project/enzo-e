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


private: // functions


private: // attributes

};

#endif /* DATA_PARTICLE_DATA_HPP */

