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
		   const EnzoConfig * enzo_config,
		   int max_initial_level) throw();

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
  virtual ~EnzoInitialMusic() throw()
  { }

  /// CHARM++ Pack / Unpack function
  void pup (PUP::er &p);
  
  /// Initialize a Block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

protected: // functions

  void throttle_input_();
  
  template <class T>
  void copy_field_data_to_array_
  (enzo_float * array, T * data,
   int mx,int my,int mz,int nx,int ny,int nz,int gx,int gy,int gz,int n4[4],
   int IX, int IY) const;

  template <class T, class S>
  void copy_particle_data_to_array_
  (T * array, S * data,
   Particle particle, int it, int ia, int np);

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

  /// Counter for throttling input
  int count_;

  /// Number of reads per process to throttle
  int throttle_count_;
  /// Number of groups with different offsets
  int throttle_offset_;
  /// Number of seconds to throttle
  int throttle_seconds_;

};

#endif /* ENZO_ENZO_INITIAL_MUSIC_HPP */

