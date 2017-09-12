// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialMusic.hpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2017-06-23
/// @brief    [\ref Enzo] Declaration of the EnzoInitialMusic class

#ifndef ENZO_ENZO_INITIAL_MUSIC_HPP
#define ENZO_ENZO_INITIAL_MUSIC_HPP

class EnzoInitialMusic : public Initial {

  /// @class    EnzoInitialMusic
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Read initial conditions from the MUSIC HDF5 files

public: // interface

  /// Constructor
  EnzoInitialMusic(int cycle,
		  double time,
		  const EnzoConfig * enzo_config) throw();

  /// Constructor
  EnzoInitialMusic() throw()
  { }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialMusic);

  /// CHARM++ migration constructor
  EnzoInitialMusic(CkMigrateMessage *m)
    : Initial (m),
      level_(0)
  {  }

  /// Destructor
  ~EnzoInitialMusic() throw()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Initialize a Block
  virtual void enforce_block
  ( Block            * block, 
    const FieldDescr * field_descr,
    const ParticleDescr * particle_descr,
    const Hierarchy  * hierarchy
    ) throw();

protected: // functions


protected: // attributes

  // NOTE: change pup() function whenever attributes change

  // Only initialize Blocks at this level
  int level_;
  
  std::vector < std::string > field_files_;
  std::vector < std::string > field_datasets_;
  std::vector < std::string > field_coords_;
  std::vector < std::string > field_names_;

  std::vector < std::string > particle_files_;
  std::vector < std::string > particle_datasets_;
  std::vector < std::string > particle_coords_;
  std::vector < std::string > particle_types_;
  std::vector < std::string > particle_attributes_;

};

#endif /* ENZO_ENZO_INITIAL_MUSIC_HPP */

