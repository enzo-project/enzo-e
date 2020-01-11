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

  /// If internode throttling enabled, sleep (i_noden * throttle_seconds_stagger_) seconds
  /// before first file open for each pe in node i_node
  void throttle_stagger_();
  /// If internode throttling enabled, sleep throttle_seconds_delay_ seconds after
  /// each open/close pair
  void throttle_delay_();
  
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

  /// Use Charm++ mutex to limit open files to one per node
  /// REQUIRES CONFIG_SMP_MODE
  bool throttle_intranode_;

  /// Throttle output between nodes by introducing a delay before
  /// starting reading based on node id throttle_group_size, and
  /// throttle_seconds_delay_
  bool throttle_internode_;
  
  /// Number of groups with different group_sizes
  int throttle_group_size_;
  /// if internode throttling, start reading only after stagger * K
  /// seconds for pe's in node K
  double throttle_seconds_stagger_;
  /// if internode throttling, delay after each open/close pair
  double throttle_seconds_delay_;
};

#endif /* ENZO_ENZO_INITIAL_MUSIC_HPP */

