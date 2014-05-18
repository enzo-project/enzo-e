// See LICENSE_CELLO file for license and copyright information

/// @file     particle_ParticleBlock.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     yyyy-mm-dd
/// @brief    [\ref Particle] Declaration of the ParticleBlock class
///

#ifndef PARTICLE_PARTICLE_BLOCK_HPP
#define PARTICLE_PARTICLE_BLOCK_HPP

class ParticleBlock {

  /// @class    ParticleBlock
  /// @ingroup  Particle
  /// @brief    [\ref Particle] 

public: // interface

  /// Constructor
  ParticleBlock();

  /// Destructor
  ~ParticleBlock();

  /// Copy constructor
  ParticleBlock(const ParticleBlock & ParticleBlock);

  /// Assignment operator
  ParticleBlock & operator= (const ParticleBlock & ParticleBlock);

  /// CHARM++ Pack / Unpack function
  inline void pup (PUP::er &p)
  {
    TRACEPUP;
    // NOTE: change this function whenever attributes change
  }
  
  void position (double * x, double * y=0, double * z=0) const {};
private: // functions


private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Positions relative to containing grid
  std::vector<short> p_[3];

  /// Scaled velocities such that new position = p_ + v_
  std::vector<short> v_[3];

};

#endif /* PARTICLE_PARTICLE_BLOCK_HPP */

