// See LICENSE_CELLO file for license and copyright information

/// @file     data_ParticleDescr.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2015-10-12

/// @brief    [\ref Data] Declaration of the ParticleDescr class

#ifndef DATA_PARTICLE_DESCR_HPP
#define DATA_PARTICLE_DESCR_HPP

class ParticleDescr {

  /// @class    ParticleDescr
  /// @ingroup  Data
  /// @brief    [\ref Data] 

  //----------------------------------------------------------------------

public: // interface

  /// Constructor
  ParticleDescr() throw();

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);

private: // functions

  //----------------------------------------------------------------------

private: // attributes

};

#endif /* DATA_PARTICLE_DESCR_HPP */

